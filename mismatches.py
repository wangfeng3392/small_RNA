#!/usr/bin/env python3
# This script calculates the mismatch rates at every nucleotide position in a given size group of small RNAs
# It requires a bam file (small RNA-seq reads alignment) and a GFF3 file (small RNA clusters)

import argparse
from itertools import islice
import pysam 
from pyfaidx import Fasta
import re

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--gff3', required=True, type=str)
	parser.add_argument('--bam', required=True, type=str)
	args = parser.parse_args()
	
	loci = hetsiRNA_loci(args.gff3)
	mismatch_counter(args.bam, loci)
	mismatch_counter_global(args.bam)


class Locus:
	def __init__(self, chrom, start, end):
		self.chr = chrom
		self.start = start
		self.end = end
	
def hetsiRNA_loci(gff):
	# Retrieve hetsiRNA (DicerCall 23 and 24 and non-miRNA) loci
	loci = {}
	with open(gff,'r') as f:
		next(f) #skip the first line in the gff3 file
		for line in f:
			info = line.strip().split('\t')
			ID = info[8].split(';')[0].split('=')[1]
			DicerCall = info[8].split(';')[1].split('=')[1]
			miRNACall = info[8].split(';')[2].split('=')[1]
			if miRNACall == "N" and DicerCall in ["23","24"]:
				chrom = info[0]
				start = int(info[3])
				end = int(info[4])
				loci[ID] = Locus(chrom, start, end)
	return loci


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

class Mismatch(object):
	def __init__(self, name, totalRead, mismatchedNT, mismatchedRead):
		self.name = name
		self.totalRead = totalRead
		self.mismatchedNT = mismatchedNT
		self.mismatchedRead = mismatchedRead

def initial_values():
	totalRead = {}
	mismatchedNT = {}
	mismatchedRead = {}
	for length in list(range(20,27)):
		totalRead[length] = 0
		mismatchedNT[length] = {}
		mismatchedRead[length] = {}
		for position in list(range(1,length+1)):
			mismatchedNT[length][position] = {"A":0, "T":0, "C":0, "G":0, "N":0}
			mismatchedRead[length][position] = 0
	return(totalRead, mismatchedNT, mismatchedRead)

def mismatch_counter(bam, loci):
	# loci is a dictionary: key=loci_ID, value=object instance in Locus class
	libraries = {}
	samfile = pysam.AlignmentFile(bam,'rb')
	for locus in loci:
		chrom = loci[locus].chr
		start = loci[locus].start - 1
		end = loci[locus].end
		for read in samfile.fetch(chrom, start, end):
			read_size = read.query_length
			read_start = int(read.reference_start) # 0-based indexing
			read_end = int(read.reference_end)
			mp = read.get_tag('MD') # mismatched position
			library = read.get_tag('RG')

			if library not in libraries:
				totalRead, mismatchedNT, mismatchedRead = initial_values()
				libraries[library] = Mismatch(library, totalRead, mismatchedNT, mismatchedRead)

			if read_size in [20, 21, 22, 23, 24, 25, 26]:
				libraries[library].totalRead[read_size] += 1
						
				interval = re.split(r'[ATGC]+', mp)
				interval = [int(x) for x in interval]
				k = len(interval)
				i = 1

				if read.is_reverse:
				# If the read is mapped to the reverse complement strand
					interval.reverse() 
					# Reverse the order of intervals
					seq = rev_comp(read.query_sequence)
				else:
					seq = read.query_sequence
				
				i = 1
				while i < k:
					readPos = sum(interval[:i]) + i
					# Transform interval lengths to absolute positions
					libraries[library].mismatchedRead[read_size][readPos] += 1
					libraries[library].mismatchedNT[read_size][readPos][seq[readPos-1]] += 1
					i += 1
					
	with open("mismatched_position.txt", 'w') as f:
		f.write(f"Library\tSize\tPosition\tRatio\tTotalCount\n")
		for library in libraries:
			for read_size in list(range(20,27)):
				for readPos in list(range(1,read_size+1)):
					ratio = libraries[library].mismatchedRead[read_size][readPos]/libraries[library].totalRead[read_size]
					f.write(f"{library}\t"
							f"{read_size}\t"
							f"{readPos}\t"
							f"{ratio:.3f}\t"
							f"{libraries[library].totalRead[read_size]}\n")

	with open("mismatched_nucleotide.txt", 'w') as f:
		f.write(f"Library\tSize\tNucleotide\tRatio\n")
		for library in libraries:
			for read_size in list(range(20,27)):
				for nucleotide in ["A", "T", "G", "C", "N"]:
					NT = libraries[library].mismatchedNT[read_size][read_size][nucleotide]
					allNT = sum(libraries[library].mismatchedNT[read_size][read_size].values())
					ratio = NT/allNT
					f.write(f"{library}\t" 
							f"{read_size}\t"
							f"{nucleotide}\t"
							f"{ratio:.3f}\n")


def mismatch_counter_global(bam):
	# loci is a dictionary: key=loci_ID, value=object instance in Locus class
	libraries = {}
	samfile = pysam.AlignmentFile(bam,'rb')
	for read in samfile:
		if not read.is_unmapped:
			read_size = read.query_length
			read_start = int(read.reference_start) # 0-based indexing
			read_end = int(read.reference_end)
			mp = read.get_tag('MD') # mismatched position
			library = read.get_tag('RG')

			if library not in libraries:
				totalRead, mismatchedNT, mismatchedRead = initial_values()
				libraries[library] = Mismatch(library, totalRead, mismatchedNT, mismatchedRead)

			if read_size in [20, 21, 22, 23, 24, 25, 26]:
				libraries[library].totalRead[read_size] += 1
						
				interval = re.split(r'[ATGC]+', mp)
				interval = [int(x) for x in interval]
				k = len(interval)
				i = 1

				if read.is_reverse:
				# If the read is mapped to the reverse complement strand
					interval.reverse() 
					# Reverse the order of intervals
					seq = rev_comp(read.query_sequence)
				else:
					seq = read.query_sequence
				i = 1
				while i < k:
					readPos = sum(interval[:i]) + i
					# Transform interval lengths to absolute positions
					libraries[library].mismatchedRead[read_size][readPos] += 1
					libraries[library].mismatchedNT[read_size][readPos][seq[readPos-1]] += 1
					i += 1
					
	with open("mismatched_position_global.txt", 'w') as f:
		f.write(f"Library\tSize\tPosition\tRatio\tTotalCount\n")
		for library in libraries:
			for read_size in list(range(20,27)):
				for readPos in list(range(1,read_size+1)):
					ratio = libraries[library].mismatchedRead[read_size][readPos]/libraries[library].totalRead[read_size]
					f.write(f"{library}\t"
							f"{read_size}\t"
							f"{readPos}\t"
							f"{ratio:.3f}\t"
							f"{libraries[library].totalRead[read_size]}\n")

	with open("mismatched_nucleotide_global.txt", 'w') as f:
		f.write(f"Library\tSize\tNucleotide\tRatio\n")
		for library in libraries:
			for read_size in list(range(20,27)):
				for nucleotide in ["A", "T", "G", "C", "N"]:
					NT = libraries[library].mismatchedNT[read_size][read_size][nucleotide]
					allNT = sum(libraries[library].mismatchedNT[read_size][read_size].values())
					ratio = NT/allNT
					f.write(f"{library}\t" 
							f"{read_size}\t"
							f"{nucleotide}\t"
							f"{ratio:.3f}\n")

if __name__ == "__main__":
	main()

