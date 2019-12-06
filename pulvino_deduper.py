#!/usr/bin/env python3

import sys
#print(version)
import argparse
import re


def get_args():
	parser = argparse.ArgumentParser(add_help=False, description="Remove PCR duplicate records on account of UMI, chromosome, soft clipping, and flag")
	parser.add_argument("-f", "--file", help="SAM input file", required=True)
	parser.add_argument("-p", "--paired", help="designates file is paired end (not single-end)", default = False, choices=['True','False'])
	parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs", default = 'STL96.txt')
	parser.add_argument("-o", "--out_prefix", help="Prefix for output files", required=False, default = '')
	parser.add_argument("-h", "--help", action = 'help', help="length of read; the index or biological read length", default=argparse.SUPPRESS)
	return parser.parse_args()
	
args = get_args()
f = args.file
o = args.out_prefix
u = args.umi
p = args.paired
if p == 'True':
	print("Not optimized for paired end yet")
	sys.exit()

UMIs = []
duplicheck = {}

chromosome_ct = 1
pcr_duplicates = 0
unique_records = 0

duplicate_trash = open(o+'duplicate_trash.sam', 'w')
deduped_apulvino = open(o+'apulvino_deduped.sam', 'w')


def ParseCigar(CIGAR, POS):
	'''Takes in a CIGAR string and returns 5' position of the read.'''
	seq = re.search("^(\d+)S", CIGAR)
	if seq == None:
		SoftClip = 0
		adjpos = POS - SoftClip
	else:
		SoftClip = int(seq.group(0)[:-1])
		adjpos = POS - SoftClip
	return adjpos
		
def ParseCigarMinus(CIGAR, POS):
	'''Takes in a CIGAR string and returns 5' position of the read for the minus strand.'''
	seq = re.findall("\d+S$|\d+N|\d+D|\d+M", CIGAR)
	#print(CIGAR, seq)
	adjpos = 0
	for i in seq:
		adjpos = POS + int(i[:-1])
	
	return adjpos


with open(u, "r") as u:
    for line in u:
        line = line.strip()
        UMIs.append(line)
        

with open(f, "r") as file:
	for line in file:
		if '@' == line[0]:
			deduped_apulvino.write(line)
		elif line[0] != '@':
			line = line.strip().split()
			CHROM = line[2]
			POS = int(line[3])
			CIGAR = line[5]
			flag = int(line[1])
			umi = re.search("[A-Z]+$", line[0]).group(0)
			if umi not in UMIs:
				#print(umi)
				#print(umi not in UMIs)
				continue
			if (flag & 16 ) != 16:
				strand = "plus"
				adjpos = ParseCigar(CIGAR, POS)
			else:
				strand = "minus"
				adjpos = ParseCigarMinus(CIGAR, POS)
			
			
			#print(umi)
			
			
			
			dupetupe = (strand, umi, adjpos, CHROM)
			#reads map to plus strand
			if dupetupe not in duplicheck:
				unique_records += 1
				duplicheck[dupetupe] = line
				deduped_apulvino.write("\t".join(line)+"\n")
				#deduped_apulvino.write(line)
				#deduped_apulvino.write('\n')
					
			else:
				pcr_duplicates += 1
				duplicate_trash.write("\t".join(line)+"\n")
				
	print("Number of PCR duplicates:", pcr_duplicates)		
	print("Number of Unique Clones:", unique_records)								
 		
				
				
				
				
deduped_apulvino.close()		
duplicate_trash.close()		
	
		