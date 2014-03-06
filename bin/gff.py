#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that from a gene id of GFF file and specified length, gives you the number to extract the upstream sequence of this gene

# Imports
from Bio import SeqIO
from BCBio import GFF
import time

def retrievePos(geneid,gfffile,fastafile,length=100):
	"""
	This function takes a gene id corresponding to a specific gene in a GFF file, and a length (default = 100nt), it retrieves the start and end position and verify the eventual overlap
	"""

	rec = [] # list of records of each entry in gff file

	# GFF file parser
	with open(gfffile,"r") as f:
		for line in GFF.parse(f):
			rec.append(line) # scrap each id in GFF file and put it in a list

	# Looping into content to identify start and end
	postable = [] # table of start, end, strand and sequence data on each genes
	for r in rec:
		feat = r.features # features of each record
		for f in feat: # looping in all features (genes) of each sequence (scaffold)
			name = f.id
			start = f.location.start.position # retrieve int as position
			end = f.location.end.position
			strand = f.location.strand

			if name == geneid:
				seqid = r.id
				genestart = start
				geneend = end
				genestrand = strand
			
			postable.append([start,end,strand,r.id]) #constructing a position matrix


	# Verifying there is no overlap between upstream sequence you to extract and anothergenes
	
	# Remove positions of sequences differents from sequences of the gene of interest
	I = set()
	for pos in postable:
		if pos[3] != seqid:
			I.add(postable.index(pos))
	for i in sorted(I,reverse=True):
		del pos[i]
			
	
	postable.sort() #sort for starting position

	upstream = []
	if genestrand == 1:
		extract = genestart - length # beginning of the extract
		if extract <= 0: # we don't want to extract negative bases!
			extract = 1
		for pos in postable:
			if pos[1] in range(extract,genestart):# if the gene is overlapping
				extract = pos[1]+1
		upstream = [extract,genestart-1,genestrand] # gives the sequence you want to extract upstream of the gene
	elif genestrand == -1:
		extract = geneend + length # end of the extract
		seqlen = fastalen(fastafile,seqid)
		if extract >= fastalen : #if extract is over the sequence
			extract = seqlen
		for pos in postable:
			if pos[0] in range(geneend+1,extract+1):
				extract = pos[0]-1
		upstream = [geneend+1,extract,genestrand]

	return upstream

def fastalen(fastafile,seqid):
	"""
	This function returns the length of fasta sequence using the sequence id and the fastafile
	"""
	length = 0

	# Find the length in fastafile of seqid
	with open(fastafile,"r") as f:
		for record in SeqIO.parse(f,"fasta"): #parse sequence after sequence the fasta file
			if record.id == seqid:
				length = len(record.seq) #length of the sequence

	if length == 0:
		print "Sequence {} was not found in {}".format(seqid,fastafile)
	return length



