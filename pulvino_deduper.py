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
#arg parse options for imported SAM file for deduplication, a paired end option which outputs error message as the program is not yet optimized
#for paired end data, -u argument so the user can pass a list of known umis used, -o argument to specify custom output file prefixes,
#and a -h argument offering a help message
UMIs = []
duplicheck = {}
#initializes UMI list and duplicheck dictionary

chromosome_ct = 1
pcr_duplicates = 0
unique_records = 0
#initialized counters for keeping track of chromosomes the program encounters, pcr_duplicate records encounterd, and unique_records encountered

duplicate_trash = open(o+'duplicate_trash.sam', 'w')
deduped_apulvino = open(o+'apulvino_deduped.sam', 'w')
#opens output files for duplicate records (trash) and for deduplicated records to keep

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
#function for parsing plus strand CIGAR strings

def ParseCigarMinus(CIGAR, POS):
	'''Takes in a CIGAR string and returns 5' position of the read for the minus strand.'''
	seq = re.findall("\d+S$|\d+N|\d+D|\d+M", CIGAR)
	#print(CIGAR, seq)
	adjpos = 0
	for i in seq:
		adjpos = POS + int(i[:-1])
	
	return adjpos
#function for parsing minus strand CIGAR strings

with open(u, "r") as u:
    for line in u:
        line = line.strip()
        UMIs.append(line)
#opens umi file passed and if umi is in the line it adds it to the UMI list initialized above

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
#opens SAM input file and if the header line is encountered these lines are written to an output file. when a record is reached
#the lines are stripped and split and variables are set for Chromosome, position, CIGAR string, flag, and Umi to be kept track of later
#if umis are not in the list return to the top of the loop until an umi from the list (one actually existing that is not unknown) is encountered
#when known umis are encountered and the plus strand is encountered call strand "plus" and invoke the plus strand CIGAR parsing function
#Else, call the strand the minus strand and invoke the minus strand CIGAR parsing function
			
			#print(umi)
			
			
			
			dupetupe = (strand, umi, adjpos, CHROM)
			#reads map to plus strand
			if dupetupe not in duplicheck:
				unique_records += 1
				duplicheck[dupetupe] = line
				deduped_apulvino.write("\t".join(line)+"\n")
				#deduped_apulvino.write(line)
				#deduped_apulvino.write('\n')
#set dupetupe, a tuple, to include information regarding strand, umi, adjusted position, and chromosome per record.
#if the tuple of this information per a certain record is not in the dictionary then count this a unique record and increment 
#unique_records, add this tuple of unique info (of strand, umi, adjusted position, and chromosome per record) to the dictionary
#and write that unique record to an output file of deduplicated output records
			else:
				pcr_duplicates += 1
				duplicate_trash.write("\t".join(line)+"\n")
#Otherwise, the record is a pcr duplicate and the counter is incremented to keep count of the given duplicate record, and that record
#is written out to a file containing all duplicate records encountered for the user's reference
	print("Number of PCR duplicates:", pcr_duplicates)		
	print("Number of Unique Clones:", unique_records)								
#the information recorded in standard out will be accompanied by a recorded number of PCR duplicates this program encountered and
#the number of unique records encountered as well for the user's reference
				
				
				
				
deduped_apulvino.close()		
duplicate_trash.close()
#outfiles are closed upon the program's completed accounting for all PCR duplicate and non-duplicate
#records from the SAM files passed to the program by the user
	
		
