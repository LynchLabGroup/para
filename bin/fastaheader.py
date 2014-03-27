#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This scripts rewrite Fasta header to have sequence name in first place

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def shorten_name(seqrecord):
	"""This function return shortened name of given sequence."""

	name = seqrecord.description

	name = name.split("_|_")

	ident = name[1]

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

def rewrite_fasta(filename,output=None):
	"""Take a fasta file and rewrite a fasta file with shortened name."""
	if output == None:
		output = filename.split(".")[0] + "rw.fasta"

	rec_out = []
	with open(filename,"r") as f:
		for record in SeqIO.parse(f,"fasta"):
			newname = shorten_name(record)

			rec = SeqRecord(record.seq,id=newname,name=newname,description="")
			rec_out.append(rec)
	
	print "Writing file {}".format(output)
	SeqIO.write(rec_out,output,"fasta")
	print "Done."