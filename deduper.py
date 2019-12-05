#!/usr/bin/env python3

from sys import version
import gzip
#print(version)
import argparse
import re


def get_args():
	parser = argparse.ArgumentParser(description="organize table of renamed file and if it's index or biological file")
	parser.add_argument("-f", "--file", help="SAM input file", required=True)
	parser.add_argument("-p", "--paired", type=int, help="designates file is paired end (not single-end)", required=False)
	parser.add_argument("-u", "--umi", type=int, help="designates file containing the list of UMIs", required=False)
	parser.add_argument("-h", "--help", type=int, help="length of read; the index or biological read length", required=False)
	return parser.parse_args()
	
args = get_args()
f = args.file
p = args.paired
u = args.umi
h = args.help

UMIs = []
duplicheck = {}

chromosome_ct = 1


duplicate_trash = open('duplicate_trash.sam', 'w')
deduped_apulvino = open('deduped_apulvino.sam', 'w')
unmapped = open('unmapped.sam', 'w')

def ParseCigar(CIGAR, POS):
	'''Takes in a CIGAR string and returns 5' position of the read.'''
	seq = re.search("^(\d+)S", CIGAR)
	if seq == None:
		SoftClip = 0
		adjpos = POS - SoftClip
	else:
		SoftClip = int(sequence.group(0))[:-1]
		adjpos = POS - SoftClip
	return adjpos
		
def ParseCigarMinus(CIGAR, POS):
	'''Takes in a CIGAR string and returns 5' position of the read for the minus strand.'''
	seq = re.search("\d+S$|\d+N|\d+D|\d+M", CIGAR)
	if seq == None:
		SoftClip = 0
		adjpos = POS - SoftClip
	else:
		SoftClip = int(sequence.group(0))[:-1]
		adjpos = POS - SoftClip
	return adjpos

fh = "/projects/bgmp/apulvino/Bi624/Deduper/STL96.txt"

f = "/projects/bgmp/apulvino/Bi624/Deduper/test.sorted.sam"

with open(fh, "r") as fh:
    for line in fh:
        line = line.strip()
        UMIs.append(line)
        

with open(f, "r") as file:
	for line in file:
		if '@' in line[0]:
			deduped_apulvino.write(line)
		elif line[0] != '@':
			line = line.strip().split()
			CHROM = line[2]
			POS = int(line[3])
			CIGAR = line[6]
			flag = int(line[1])
			umi = re.search("[A-Z]+$", line[0]).group(0)
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
				duplicheck[dupetupe] = line
				deduped_apulvino.write("	".join(line)+"\n")
				#deduped_apulvino.write(line)
				#deduped_apulvino.write('\n')
					
			else:
				duplicate_trash.write("    ".join(line)+"\n")
				#duplicate_trash.write(line)
				#duplicate_trash.write('\n')
				#if chromosome_ct != line[2]:
				#	chromosome_ct = line[2]
				#	duplicheck = {}
					
						
				

			

 		#if line[0] != '@':
 		#	line = line.strip().split()
 		#	CHROM = line[2]
 		#	POS = int(line[3])
 		#	CIGAR = line[6]
 		#	flag = int(line[1])
 		#	umi = re.search("[A-Z]+$", line[0]).group(0)


			
			#adjpos = ParseCigarMinus(CIGAR, POS)
			
			#dupetupe = (((flag & 16 ) == 16), umi, adjpos, CHROM)
			#duplicheck[dupetupe] = line
					
			#if adjpos and CHROM and umi and ((flag & 16 ) == 16) not in duplicheck:
			#	duplicheck[dupetupe] = line
			#else:
			#	deduped_apulvino.write("	".join(line)+"\n")
		#	if chromosome_ct != line[2]:
		#		chromosome_ct = line[2]
		#		for dupetupe in duplicheck:
		#			deduped_apulvino.write("    ".join(line)+"\n")
		#			duplicheck = {}
			#	if dupetupe in line:
			#		deduped_apulvino.write(line)
			#		deduped_apulvino.write('\n')
							
		#	else:
		#		duplicate_trash.write("    ".join(line)+"\n")
		#		duplicheck = {}
		#	if dupetupe in line:
		#		duplicate_trash.write(line)
		#		duplicate_trash.write('\n')
				
			#elif ((flag & 4) != 4): 
				#an unmapped read!
				#unmapped.write(line)
				#unmapped.write('\n')
				
				
				
				
deduped_apulvino.close()		
duplicate_trash.close()		
unmapped.close()	
		