#!/usr/bin/env python3

import argparse
from itertools import islice
import pysam 
from pyfaidx import Fasta
import re
import os
import subprocess

def main():
	parser = argparse.ArgumentParser()

	parser.add_argument('--size', action="store_true", default=False, \
		help="For calculating general size distribution of RNAs")

	parser.add_argument('--bam', type=str, help="Alignment file, bam format")
	parser.add_argument('--gff', type=str, help="small RNA cluster GFF file")

	options = parser.parse_args()
	if not options.size:
		print("Error: no actions specified. Try '--size' to calculate the normalized sRNA counts")
	else:
		if not options.bam:
			print("Error: no bam file specified")
		elif not options.gff:
			print("Error: no gff file of sRNA clustering specified")
		else:
			total_count = mapped_read_counter(options.bam)
			loci = sRNA_loci(options.gff)
			sRNA_size_distribution(total_count, loci, options.bam)
	
def mapped_read_counter(bamfile):
	# Mapped read count for each library in a BAM file
	# A BAM file can contain reads from different libraries
	total_count = {}
	samfile = pysam.AlignmentFile(bamfile,'rb')
	for read in samfile:
		if read.is_unmapped == False:
			library = read.get_tag('RG')
			if library in total_count:
				total_count[library] += 1
			else:
				total_count[library] = 1
	return total_count

def sRNA_loci(gff):
	# Retrieve hetsiRNA (dicerCall 23 and 24 and non-miRNA) loci
	loci = []
	with open(gff,'r') as f:
		next(f) #skip the first line in the gff3 file
		for line in f:
			info = line.strip().split('\t')
			ID = info[8].split(';')[0].split('=')[1]
			#dicerCall = info[8].split(';')[1].split('=')[1]
			#miRNACall = info[8].split(';')[2].split('=')[1]
			locus = info[0] + ":" + info[3] + "-" +info[4]
			#if miRNACall == "N" and dicerCall in ["23","24"]:
			loci.append(locus)
	return loci

def sRNA_size_distribution(total_count, loci, bamfile):
	# prepare a dictionary for siRNA counts
	siRNA_count = {}
	for library in total_count:
		siRNA_count[library] = {}

	# test if the index file (.bai) is in the directory
	bam_index = bamfile + ".bai"
	if not os.path.isfile(bam_index):
		print("No index file found. Preparing a index file for the BAM file")
		subprocess.call("samtools index " + bamfile, shell=True)
	
	# start searching sRNAs in hetsiRNA loci
	samfile = pysam.AlignmentFile(bamfile,'rb')
	for locus in loci:
		chrom = locus.split(':')[0]
		start = int(locus.split(':')[1].split('-')[0]) - 1
		end = int(locus.split(':')[1].split('-')[1])
		for read in samfile.fetch(chrom,start,end):
			# loop through all alignments in this locus
			read_size = int(read.query_length)
			library = read.get_tag('RG')
			if read_size not in siRNA_count[library]:
				siRNA_count[library][read_size] = 1
			else:
				siRNA_count[library][read_size] += 1

	# start normalization and output to a file
	with open('normalized_sRNA_distribution.txt','w') as f:
		f.write(f"Library\tSize\tRaw_count\tTotal_count\tNormalized_count\n")
		for library in siRNA_count:
			for read_size in sorted(siRNA_count[library].keys()):
				cpm = siRNA_count[library][read_size]/total_count[library]*1000000
				f.write(f"{library}\t")
				f.write(f"{read_size}\t")
				f.write(f"{siRNA_count[library][read_size]}\t")
				f.write(f"{total_count[library]}\t")
				f.write(f"{cpm:.4f}\n")


if __name__ == "__main__":
	main()








