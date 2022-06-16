#!/usr/bin/env python
# This script can simulate het-siRNA biogenesis based provided gff3 file
# In general, it takes about 1.5 seconds to generate data for 1 locus.

import sys
import subprocess
import os
import random 
from time import time
import datetime
import re
import argparse
import multiprocessing as mp
import pysam
from pyfaidx import Fasta
import numpy as np
import logging


def main():

	start_time = time()
	parser = argparse.ArgumentParser()
	parser.add_argument('--genome', type=str, help="Reference genome file, fasta format")
	parser.add_argument('--bam', type=str, help="Alignment file, bam format")
	parser.add_argument('--gff3', type=str, help="small RNA cluster GFF file")
	parser.add_argument('--output', type=str, help="prefix of output files")
	parser.add_argument('--path', action="store_true", default=False)
	
	opt = parser.parse_args()

	# prepare a log file
	#logfile = opt.output + ".log"
	#logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",
	#	format="%(asctime)-15s %(levelname)-8s %(message)s")

	# check dependencies
	dependencies(opt)

	# preparing loci
	loci = het_siRNA_loci(opt.gff3)
	
	# analyzing small RNA size distribution in all loci
	bins = size_distr(opt.bam, loci)
	
	# simulating:
	for locus in bins:
		read_key_list = get_read_keys(bins, locus)
		simulate_seqs(opt.genome, opt.bam, read_key_list, opt.output)

	end_time = time()
	
	print("Time used %.2f seconds" % (end_time - start_time))

def dependencies(opt):

	# this function validates the existance of a path
	if opt.path:
		if not os.path.exists(opt.path):
			print ("Exit: Invalid path provided.\n")
			sys.exit()
	else:
		print ("Output files will be write to working directory.\n")

	print ("Checking dependencies:")

	FNULL = open(os.devnull, 'w')
	# is bowtie installed?
	try:
		subprocess.check_call("which bowtie", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	except:
		print ("Exit: bowtie v1 is required to run this program.")
		sys.exit()

	# is samtools installed?
	try:
		subprocess.check_call("which samtools", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	except:
		print ("Exit: samtools is required to run this program.")
		sys.exit()

	# is there a file named 'simulated_smallRNAs.fa' already?
	out_fa = opt.output + ".fa"
	if os.path.exists(out_fa):
		print ('simulated_smallRNAs.fa exits already!')
		sys.exit()

	# is there a file named 'simulation_origins.sam' already?
	out_sam = opt.output + ".sam"
	if os.path.exists(out_sam):
		print ('simulation_origins.sam exits already!')
		sys.exit()

	# Does the genome file exist?
	if os.path.exists(opt.genome):
		print (f"{opt.genome} --> Found")
	else:
		print (f"{opt.genome} --> Exit: genome not found")
		sys.exit()

	# Is the genome indexed by Bowtie 1?
	bowtie_index = opt.genome.split(".")[0] + ".1.ebwt"
	if os.path.exists(bowtie_index):
		print ("Reference genome is indexed")
	else:
		print ("Reference genome is not indexed")
		print ("start indexing with 'bowtie-build'...")

		# indexing the genome with bowtie-build
		subprocess.call(["bowtie-build", opt.genome, opt.genome.split(".")[0]], 
			shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

		if os.path.exists(bowtie_index):
			print ("Reference genome indexed!")
		else:
			print ("FAILED to index the genome with bowtie-build!")
			sys.exit()

	# Is there a .fai file?
	genome_fai = opt.genome + ".fai"
	if os.path.exists(genome_fai):
		print ("Genome .fai index found")
	else:
		print ("Genome index file (.fai) not found")
		time.sleep(1)
		print ("Start building .fai file with samtools...")

		# making a fai file with samtools faidx
		subprocess.call(["samtools faidx", genome], shell=True)

		if os.path.exist(genome_fai):
			print ("Reference fai file created!")
		else:
			print ("FAILED to create a fai file with samtools!")
			sys.exit()

	# test if the index file (.bai) is in the directory
	
	bam_index = opt.bam + ".bai"
	if not os.path.exists(bam_index):
		print ("No index file found. Preparing a index file for the BAM file")
		subprocess.call("samtools index " + bam, shell=True, 
			stdout=FNULL, stderr=subprocess.STDOUT)
		print ("Sample BAM file indexed!")
	else:
		print("Sample BAM file indexed!")
	print ("\n")

def progress(total_reads, read_count):
	# this function is called to update progress every 1000 reads
	if read_count % 1000 == 0:
		p = round(read_count/total_reads*20)
		q = read_count/total_reads*100
		sys.stdout.write('\r')
		sys.stdout.write("[%-20s]%d%% " % ('='*p, q))
		sys.stdout.flush()
	else:
		pass

def get_new_base(old_base):
	# select a new base for the simulated mismatched nucleotide in a given read
	# choose a base other than itself

	nucleotides = ["A", "C", "T", "G"]
	try:
		nucleotides.remove(old_base.upper())
		new_base = random.choice(nucleotides)
	except ValueError:
		# if old_base is not in A, C, T and G
		new_base = "N"
	return new_base

def seq_err(seq):
	# simulate sequencing errors: incorporate mismatches into sequences
	# options: 1 nt or 2 nt mismatches allowed in a read
	# default per-read probability of 1 error: 0.15
	# default per-read probability of 2 errors: 0.05
	# probability of mismatch(es) are high because of variations between ecotypes
	
	pick = random.randint(1, 100)
	if pick <= 15:
		# read with 1 mismatch

		mutated = random.randint(0, len(seq) - 1)
		old_base = seq[mutated]
		new_base = get_new_base(old_base)
		after = len(seq) - mutated - 1
		MD = f"MD:Z:{mutated}{new_base}{after}"

		s = list(seq)
		s[mutated] = new_base.lower()
		err_read = ''.join(s)
		err_read = (MD, err_read)

	elif pick >15 and pick <=20:
		# read with 2 mismatches
		while True:
			mutated_1 = random.randint(0, len(seq) - 1)
			mutated_2 = random.randint(0, len(seq) - 1)
			if mutated_1 != mutated_2:
				mutated_1,mutated_2 = min(mutated_1,mutated_2), max(mutated_1,mutated_2)
				break

		old_base_1 = seq[mutated_1]
		new_base_1 = get_new_base(old_base_1)

		old_base_2 = seq[mutated_2]
		new_base_2 = get_new_base(old_base_2)

		between = mutated_2 - mutated_1 - 1
		after = len(seq) - mutated_2 - 1

		MD = f"MD:Z:{mutated_1}{new_base_1}{between}{new_base_2}{after}"

		s = list(seq)
		s[mutated_1] = new_base_1.lower()
		s[mutated_2] = new_base_2.lower()
		err_read = ''.join(s)
		err_read = (MD, err_read)

	else:
		err_read = ("MD:Z:0", seq)

	return err_read

def rev_comp(seq):
	# generate reverse-complement sequence of the input

	temp = []
	comp_nt = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N",
				"a":"t", "t":"a", "g":"c", "c":"g", "n":"n"}
	for nucleotide in seq[::-1]:
		try:
			temp.append(comp_nt[nucleotide])
		except KeyError:
			temp.append("N")
	rev_comp_seq = ''.join(temp)
	return rev_comp_seq

def het_siRNA_loci(gff3):
	# this function makes a dictionary of het-siRNA loci based on the provided annotation
	# the simulated siRNA loci have the same size distribution as the provided annotation
	# genome locations of the simulated loci are randomly picked
	# overlapping of loci is prohibited

	print (f"Determining het-siRNA loci based on {gff3}\n")
	
	loci = []
	
	with open(gff3,'r') as f:
		next(f) #skip the first line in the gff3 file
		for line in f:
			info = line.strip().split('\t')
			ID = info[8].split(';')[0].split('=')[1]
			dicerCall = info[8].split(';')[1].split('=')[1]
			miRNACall = info[8].split(';')[2].split('=')[1]
			locus = info[0] + ":" + info[3] + "-" +info[4]
			if miRNACall == "N" and dicerCall in ["23","24"]:
				loci.append(locus)

	print (f"Pseudo het-siRNA loci to be prepared: {len(loci)}\n")
	
	return loci

def size_distr(bam, loci):
	# this function calculates the percentage of reads with different sizes
	# prepare empty dictionary to store sRNA size information, key represents locus location
	# value: list of 6 integers, representing the counts of 12 mers, 23 mers, 24 mers, 25 mers, 26 mers and other sizes
	
	bins = {}
	total_count = 0

	print ("Calculating counts of reads with different sizes...")

	# start searching sRNAs from sample BAM file
	samfile = pysam.AlignmentFile(bam,'rb')
	for locus in loci:
		chrom = locus.split(':')[0]
		start = int(locus.split(':')[1].split('-')[0]) - 1
		end = int(locus.split(':')[1].split('-')[1])

		# key: locus; format: Chr1:8-219
		# value: list of 6 integers -- counts of 12, 23, 24, 25, 26, other sizes
		bins[locus] = [0,0,0,0,0,0]
		
		for read in samfile.fetch(chrom,start,end):
			# loop through all alignments in this locus
			read_size = int(read.query_length)
			total_count += 1

			if read_size == 12:
				bins[locus][0] += 1
			elif read_size == 23:
				bins[locus][1] += 1
			elif read_size == 24:
				bins[locus][2] += 1
			elif read_size == 25:
				bins[locus][3] += 1
			elif read_size == 26:
				bins[locus][4] += 1
			else:
				bins[locus][5] += 1

	print ("Counting finished!")
	print (f"Pseudo het-siRNA reads to be prepared: {total_count}\n")
	
	return bins

def accumu(list):
	# calculate cummulative sum of a list
	total = 0
	for i in list:
		total += i
		yield total

def get_read_keys(bins, locus):
	# randomly generating reads for a het-siRNA locus created from sample BAM file
	# the counts of different sizes of sRNAs reflect the numbers found in the sample BAM
	# strandedness and locations of reads are randomly picked

	# n = read counts in the given locus based on sample BAM file
	read_key_list = []

	chrom = locus.split(':')[0]
	start = int(locus.split(':')[1].split('-')[0]) # 1-based offset
	end = int(locus.split(':')[1].split('-')[1]) # 1-based offset
	locus_size = end - start + 1

	# total count of sRNAs in a given locus
	total_count = sum(bins[locus])
	cummulative_sum = list(accumu(bins[locus]))

	random_array = np.random.randint(1, total_count + 1, size=(total_count,1))
	for i in random_array:
		# generate information for simulated reads
		# pick a size of sRNA read: 12, 23, 24, 25, 26, other (8-11, 13-22, 27-75)
		# np.random is much faster than random
		
		size_pick = i[0]
		if size_pick < cummulative_sum[0]:
			size = 12
		elif size_pick >= cummulative_sum[0] and size_pick < cummulative_sum[1]:
			size = 23
		elif size_pick >= cummulative_sum[1] and size_pick < cummulative_sum[2]:
			size = 24
		elif size_pick >= cummulative_sum[2] and size_pick < cummulative_sum[3]:
			size = 25
		elif size_pick >= cummulative_sum[3] and size_pick < cummulative_sum[4]:
			size = 26
		else:
			if locus_size <=27:
				size = random.choice([random.randint(8,11),random.randint(13,22)])
			elif locus_size > 27 and locus_size < 75:
				size = random.choice([random.randint(8,11),random.randint(13,22),random.randint(27,locus_size)])
			else:
				size = random.choice([random.randint(8,11),random.randint(13,22),random.randint(27,75)])

		# pick a strand where the read is from
		strand_pick = random.choice(["top","bottom"])

		# pick the left coordinate where the alignment starts
		while True:
			start_pick = random.randint(start, end)

			# read-end: start_pick + size - 1
			if start_pick + size - 1 <= end:
				break

		read_key = (chrom, strand_pick, start_pick, start_pick + size - 1)
		read_key_list.append(read_key)

	return read_key_list

def simulate_seqs(genome, bam, read_key_list, output):
	# assign sequences to simulated reads, based on read_key
	# read_key: (chrom, strand, start, end); start and end are 1-based offset
	# this function writes info of simulated read to a fasta file and a SAM file
	
	i = 1
	out_fa = output + ".fa"
	out_sam = output + ".sam"
	with open(out_fa,'a') as f1,open(out_sam,'a') as f2:
		for read_key in read_key_list:
			chrom, strand, start, end = read_key[0], read_key[1], read_key[2], read_key[3]
			region="%s:%d-%d" % (chrom, start, end) 

			# retrieve sequence by samtools: samtools faidx Athaliana_167_TAIR10.fa Chr1:1-10
			# output format: b'>Chr1:1-10\nCCCTAAACCC
			#p1= subprocess.Popen(["samtools", "faidx", opt.genome, region], stdout = subprocess.PIPE)
			#seq = p1.communicate()[0].decode().split('\n')[1:]
			#seq = ''.join(seq)

			# use pyfaidx to extract sequences from genome
			# pyfaidx uses 0-based, left inclusive, righ exlcusive offsets;
			# need to do a little fix-up for the ranges
			# pyfaidx is more than 100x faster than samtools faidx

			index_genome = Fasta(genome)
			seq = index_genome[chrom][start-1:end].seq

			# incorporate mismatch(es) by calling seq_err(), which return MD and mutated seq
			MD, seq = seq_err(seq)

			# number of mismatches: stratum
			stratum = len(re.findall(r"[ATGCN]", MD.split(":")[-1]))
			
			XA = "%s%s"% ("XA:i:", stratum)

			# reverse complement seq if read is from bottom strand
			if strand == "bottom":
				fasta_seq = rev_comp(seq)
			else:
				fasta_seq = seq

			# write info to 'simulated_smallRNAs.fa'
			# header format: >Chr1:1-20_top_MD:Z:0_1
			f1.write(f">{chrom}:{start}-{end}_{strand}_{MD}_{i}\n")
			f1.write(f"{fasta_seq}\n")

			# write info to 'smallRNAs_origins.sam'
			# name format: Chr1:1-20_top_MD:Z:0_1
			f2.write(f"{chrom}:{start}-{end}_{strand}_{MD}_{i}\t")
			f2.write(f"{0 if strand=='top' else 16}\t")
			f2.write(f"{chrom}\t")
			f2.write(f"{start}\t")
			f2.write(f"255\t")
			f2.write(f"{end - start + 1}M\t")
			f2.write(f"*\t")
			f2.write(f"0\t")
			f2.write(f"0\t")
			f2.write(f"{seq}\t")
			f2.write(f"{'I'*len(seq)}\t")
			f2.write(f"*\t")
			f2.write(f"{XA}\t")
			f2.write(f"{MD}\t")
			f2.write(f"XX:S:Simulated\t")
			f2.write(f"RG:Z:{bam}\n")

			i += 1


if __name__ == "__main__":
	main()





