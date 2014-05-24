#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This scripts rewrite Fasta header to have sequence name in first place
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def shorten_name(seqrecord,char):
	"""This function return shortened name of given sequence."""

	name = seqrecord.description

	name = name.split(char)

	ident = name[0]

	#if len(ident) > 10:
	#	ident = ident[0:10]

	return ident

def rewrite_phylip(filename,output=None):
	"""Take a fasta file and rewrite a phylip file with shortened name."""
	if output == None:
		output = filename.split(".")[0] + "rw.fasta"

	rec_out = []
	with open(filename,"r") as f:
		for record in SeqIO.parse(f,"fasta"):
			newname = shorten_name(record)

			rec = SeqRecord(record.seq,id=newname,name=newname,description="")
			rec_out.append(rec)
	
	print "Writing file {}".format(output)
	SeqIO.write(rec_out,output,"phylip-relaxed")
	print "Done."

def rewrite_fasta(filename,output=None,char=None):
	"""Take a fasta file and rewrite a fasta file with shortened name."""
	if output == None:
		output = filename.split(".")[0] + "rw.fasta"
	if char == None:
		char = "|"

	rec_out = []
	with open(filename,"r") as f:
		for record in SeqIO.parse(f,"fasta"):
			newname = shorten_name(record,char)

			rec = SeqRecord(record.seq,id=newname,name=newname,description="")
			rec_out.append(rec)
	
	print "Writing file {}".format(output)
	SeqIO.write(rec_out,output,"fasta")
	print "Done."

if __name__ == "__main__":
	### Command-line Parser ###

	parser = argparse.ArgumentParser(description="Program to simplify sequence name in fasta file.")

	parser.add_argument("fasta", help=".fasta file")

	parser.add_argument("char", help="char that separate fasta header")

	parser.add_argument("out", help="output file fasta format")

	args = parser.parse_args()

	rewrite_fasta(args.fasta,args.out,args.char)
