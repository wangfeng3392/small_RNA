#!/usr/bin/env python
# This script analyzes the spatial relationship between different sizes of small RNAs
# Only 24 nt siRNAs containing 5' A are computational reconstituted in predicted pairs;
# The predicted siRNA pairs should come from siRNA loci
# It requires a bam file (small RNA alignment) and a GFF3 file of small RNA clusters
# You can define query size of small RNAs (what size to basepair with a 5'A 24 nt siRNA)

import pysam
import sys
import os
import argparse
from time import time
import logging
import subprocess
import re

def main():
	start_time = time()
	parser = argparse.ArgumentParser()
	parser.add_argument('--bam', nargs='+', required=True, type=str, 
		help="BAM files to be analyzed")
	parser.add_argument('--gff3', type=str, help="small RNA cluster GFF file")
	parser.add_argument('--query', required=True, type=int, 
		help="size (nt) of query group: 12, 23 or 24")
	parser.add_argument('--standard', required=True, type=str)
	parser.add_argument('--output', required=True, type=str,
		help="output file name")

	opt = parser.parse_args()

	# check if all arguments are correctly provided:
	validate_args(opt)

	# analyze siRNA loci
	loci = het_siRNA_loci(opt.gff3)

	# Generate lists of 24 nt siRNAs, 23 nt siRNAs, 12 nt siRNAs
	print ("Collecting query and subject siRNA information...")
	query_collector, subject_collector = extract_siRNAs(opt, loci)

	print ("Collecting done! Start analyzing spatial relationships...")
	p5_offset(opt, query_collector, subject_collector)
	
	time_elapsed = time() - start_time
	print (f"This analysis takes {time_elapsed} seconds")


def validate_args(opt):
	# Check bam and bai files
	for bam in opt.bam:
		bam_index = bam + ".bai"

		# Are the bam files in the working directory?
		if os.path.exists(bam):
			print (f"{bam} ----> Found!")
		else:
			print (f"FATAL: {bam} is not in the working directory!")
			sys.exit()

		# Are the bai files in the working directory?
		if os.path.exists(bam_index):
			print (f"{bam_index} ----> Found!")
		else:
			print (f"{bam_index} is not in the working directory")
			print (f"Starting index {bam} with samtools...")
			
			subprocess.call(["samtools", "index", bam], shell=True)
			
			if os.path.exists(bam_index):
				print (f"{bam_index} ----> Created!")
			else:
				print (f"FATAL: failed to create the index file for {bam}, quiting now!")
				sys.exit()


	# Check argument for groups: 23vs24, 12vs24 and 24vs24
	if opt.query in [12, 23, 24]:
		print (f"This script will analyze the spatial relationship between {opt.query} nt siRNAs and 24 nt siRNAs")
	else:
		print ("Choose a size group for analysis: 12, 23 or 24")
		print ("Try: python spatial_relationship.py --bam YOUR_BAM.bam --query 12")
		sys.exit()


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

def mismatch_positions(read):
	# This function read the MD field of an alignment in a BAM/SAM file and returns a list of mismatched positions
	# format of md: 22A0

	mismatch_positions_list = []
	mismatch_nts_list = []

	md = read.get_tag("MD")

	interval = re.split(r'[A-Z]', md)
	interval = [int(x) for x in interval]

	if read.is_reverse:
		# If the read is mapped to the reverse complement strand 
		# Reverse the order of intervals, generate reverse-complement sequences
		interval.reverse()
		seq = rev_comp(read.query_sequence)
	else:
		seq = read.query_sequence
	
	i = 1
	k = len(interval)
	while i < k:
		# Transform interval lengths to absolute positions
		mpos = sum(interval[:i]) + i
		mismatch_positions_list.append(mpos)
		mismatch_nts_list.append(seq[mpos - 1])
		i += 1
	
	return (mismatch_positions_list, mismatch_nts_list)

def collecting_keys(read):
	# get key information from BAM/SAM file

	chrom = read.reference_name
	strand = '-' if read.is_reverse else '+'
	left = read.reference_start # 0-based offset
	right = read.reference_end # 0-based offset, points to one past the last aligned residue

	# 5' and 3' ends of the read
	if read.is_reverse:
		seq = rev_comp(read.query_sequence)
	else:
		seq = read.query_sequence
	p5_nt = seq[0]
	p3_nt = seq[-1]

	# mismatch positions and mismatch nt in the read
	mismatch_positions_list, mismatch_nts_list = mismatch_positions(read)

	# if the 3' last nt is a mismatch (positions are 0-based):
	try:
		if mismatch_positions_list[-1] == int(read.query_length) -1:
			p3_mismatch = mismatch_nts_list[-1]
		else:
			p3_mismatch = "0"
	except IndexError:
		p3_mismatch = "0"

	# fill in query_collector now
	collected_key = (chrom, strand, left, right, p5_nt, p3_nt, p3_mismatch)
	#collected_key = (chrom, strand, left, right)
	return collected_key

def het_siRNA_loci(gff3):
	# this function makes a dictionary of het-siRNA loci based on the provided annotation
	# the simulated siRNA loci have the same size distribution as the provided annotation
	# genome locations of the simulated loci are randomly picked
	# overlapping of loci is prohibited

	print (f"Determining het-siRNA loci based on {gff3}")
	
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

	print (f"Number of het-siRNA loci: {len(loci)}")
	
	return loci

def extract_siRNAs(opt, loci):
	# prepare empty dictionaries to collect siRNA information; key= (bam file name)
	query_collector = {}
	subject_collector = {}

	bams = opt.bam
	query_size = int(opt.query)
	for bam in bams:
		query_collector[bam] = {}
		subject_collector[bam] = {}

		samfile = pysam.AlignmentFile(bam, 'rb')
		for locus in loci:
			# Pysam is 0-based; gff3 is 1-based;
			# pysam uses left-inclusive and right exclusive: [)
			# gff3 uses left inclusive and right inclusive: []
			# start - 1 to convert
			chrom = locus.split(':')[0]
			start = int(locus.split(':')[1].split('-')[0]) - 1
			end = int(locus.split(':')[1].split('-')[1])

			for read in samfile.fetch(chrom,start,end):
				if not read.is_unmapped:
					# get basic information from SAM file

					read_size = int(read.query_length)
					if query_size != 24:
						# if query_size is 12 or 23, record first 4 values of key tuples
						if read_size == query_size:
							# if the length of this read equals to query size (12, 23):
							collected_key = collecting_keys(read)
							collected_key = collected_key[0:4]
							
							if collected_key not in query_collector[bam]:
								query_collector[bam][collected_key] = 1
							else:
								query_collector[bam][collected_key] += 1
						elif read_size == 24:
							# if read size is 24, this is subject RNA, record information in subject collector
							collected_key = collecting_keys(read)
							
							if collected_key[4] == "A":
								# if the 5' nucleotide is A, put it into subject collector
								collected_key = collected_key[0:4]
								if collected_key not in subject_collector[bam]:
									subject_collector[bam][collected_key] = 1
								else:
									subject_collector[bam][collected_key] += 1

					else:
						collected_key = collecting_keys(read)
						if collected_key[4] == "A":
							# if the 5' nucleotide is A, put it into subject collector
							collected_key = collected_key[0:4]
							if collected_key not in subject_collector[bam]:
								subject_collector[bam][collected_key] = 1
							else:
								subject_collector[bam][collected_key] += 1
						else:
							# else, put it into query collector
							collected_key = collected_key[0:4]
							if collected_key not in query_collector[bam]:
								query_collector[bam][collected_key] = 1
							else:
								query_collector[bam][collected_key] += 1


	return (query_collector, subject_collector)


def update_collectors(old_collector):
	new_collector = {}
	for bam in old_collector:
		new_collector[bam] = {}
		for key in old_collector[bam]:
			new_key = key[0:4]
			if new_key not in new_collector[bam]:
				new_collector[bam][new_key] = old_collector[bam][key]
			else:
				new_collector[bam][new_key] += old_collector[bam][key]

	return new_collector


def p5_offset(opt, query_collector, subject_collector):
	# prepare a dictionary to store offset values
	offset_dict = {}
	query_size = int(opt.query)

	for bam in opt.bam:
		offset_dict[bam] = {}
		for offset in range(1, query_size + 24):
			offset_dict[bam][offset] = 0
	"""
	if opt.standard == "stringent":
		# stringent searching criteria:
		# subject 24 nt siRNAs should have 5'As
		# query 23 nt siRNAs should have 3' end mismatches
		# queries and subjects are on different strands
		pass
		
	elif opt.standard == "general":
		# update collectors based on the first 4 values in the key tuples
		new_query_collector = update_collectors(query_collector)
		new_subject_collector = update_collectors(subject_collector)
	"""
		# only requires opposite strandedness for queries and subjects
	for bam in query_collector:
		for q_key in query_collector[bam]:
			q_chrom, q_strand, q_left, q_right = q_key
			s_chrom = q_chrom
			
			for offset in range(1, query_size + 24):
				if q_strand == '-':
					s_strand = '+'
					s_left = q_right - offset
					s_right = s_left + 24
				else:
					s_strand = '-'
					s_right = q_left + offset
					s_left = s_right - 24

				s_key = (s_chrom, s_strand, s_left, s_right)
				if s_key in subject_collector[bam]:
					if offset in offset_dict[bam]:
						#offset_dict[bam][offset] += subject_collector[bam][s_key]*query_collector[bam][q_key]
						offset_dict[bam][offset] += 1
					else:
						#offset_dict[bam][offset] = subject_collector[bam][s_key]*query_collector[bam][q_key]
						offset_dict[bam][offset] = 1
	
	with open(opt.output,'w') as f:
			f.write(f"Library\tOffset\tOccurence\tFrequency\n")
			for bam in offset_dict:
				for offset in range(1, query_size + 24):
					total_occurence = sum(offset_dict[bam].values())
					perc = offset_dict[bam][offset]*100/total_occurence
					f.write(f"{bam}\t")
					f.write(f"{offset}\t")
					f.write(f"{total_occurence}\t")
					f.write(f"{perc}\n")

if __name__ == "__main__":
	main()












