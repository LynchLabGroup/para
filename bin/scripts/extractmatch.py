#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def name_set(fastafile):
	"""
	This function returns the set of name of a fasta file.
	"""
	names = set()
	with open(fastafile,"r") as f:
		for record in SeqIO.parse(f,"fasta"):
			names.add(record.id)

	return names

def select_write_fasta(names,fastafile,outputfile):
	"""
	Write a given fastafile as outputfile file using only named sequences.
	"""

	out = []

	with open(fastafile,"r") as f:
		for record in SeqIO.parse(f,"fasta"):
			if record.id in names:
				out.append(record)

	print "Writing file {}...".format(outputfile)
	SeqIO.write(out,outputfile,"fasta")
	print "Done."

def main():
	
	parser = argparse.ArgumentParser(description="Simple script to rewrite a fasta file using the names of another fasta file.")

	parser.add_argument("name", help="Fasta file from which you want to extract names")
	parser.add_argument("inp",help="Fasta file you want to rewrite")
	parser.add_argument("out",help="Name of the output fasta file")	

	args = parser.parse_args()

	names = name_set(args.name)

	select_write_fasta(names,args.inp,args.out)

if __name__ == "__main__":
	main()