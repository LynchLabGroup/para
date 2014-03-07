#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that from a gene id of GFF file and specified length, gives you the number to extract the upstream sequence of this gene

# Imports
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

def retrievePos(geneids,gfffile,fastafile,length=100):
	"""
	DESCRIPTION:

	Take a list of gene ids, and retrieve the upstream sequence (default = 100nt) from fastafile, verifying in gfffile that there is no overlap. Return a list of genes, locations and sequences.

	USAGE:

	geneids -- list of gene ids in gfffile
	gfffile -- a GFF file
	fastafile -- a Fasta file
	length (optional) -- length of the upstream sequence to extract
	"""
	rec = [] # list of records of each entry in gff file

	# GFF file parser, puts GFF file in memory
	with open(gfffile,"r") as f:
		print "Parsing file {}...".format(gfffile)
		for line in GFF.parse(f):
			rec.append(line) # scrap each id in GFF file and put it in a list

	interest = [] # list of coordinates of gene of interest

	# Looping into content to identify characteristics of gene of interest
	postable = [] # table of start, end, strand and sequence data on each genes
	print "Extracting genes of interest..."
	for r in rec:
		feat = r.features # features of each record
		for f in feat: # looping in all features (genes) of each sequence (scaffold)
			name = f.id
			start = f.sub_features[0].location.start.position # retrieve int as position
			end = f.sub_features[0].location.end.position 
			strand = f.sub_features[0].location.strand

			# Retrieves data from genes of interest
			for g in geneids:
				if name == g:
					interest.append([start,end,strand,g,r.id])
			
			postable.append([start,end,strand,r.id]) #constructing a position matrix

	# Retrieves start and end position to extract, while verifying there is no overlapping genes
	upstream = []
	print "Retrieving upstream locations..."
	for i in interest:
		start = i[0] # natural beginning of the gene (begins with position 1)
		end = i[1] # natural ending
		strand = i[2] # strand 1: plus -1: minus
		gene = i[3] # gene id
		seqid = i[4] # sequence id on which gene is located
		extract = 0
		if strand == 1:
			extract = start - length # beginning of the extract
			if extract <= 0: # we don't want to extract negative bases!
				extract = 1
			for pos in postable:
				if pos[-1] == seqid and pos[1] in range(extract,start):# if the gene is overlapping on the same scaffold
					extract = pos[1]+1
			upstream.append([extract,start,strand,gene,seqid])# gives the sequence you want to extract upstream of the gene
		elif strand == -1:
			extract = end + length # end of the extract
			seqlen = fastalen(fastafile,seqid)
			if extract >= fastalen : #if extract is over the sequence
				extract = seqlen
			for pos in postable:
				if pos[-1] == seqid and pos[0] in range(end+1,extract+1):
					extract = pos[0]
			upstream.append([end+1,extract,strand,gene,seqid])

	print "Retrieving sequences..."
	upstream = retrieveSeq(fastafile,upstream)
	
	return upstream

def fastalen(fastafile,seqid):
	"""Return length of sequence seqid in fastafile."""

	length = 0

	# Find the length in fastafile of seqid
	with open(fastafile,"r") as f:
		for record in SeqIO.parse(f,"fasta"): #parse sequence after sequence the fasta file
			if record.id == seqid:
				length = len(record.seq) #length of the sequence

	if length == 0:
		print "Sequence {} was not found in {}".format(seqid,fastafile)
	return length


def retrieveSeq(fastafile,uplist):
	"""Take retrievePos return and return the same list with sequences appended."""

	records = [] # list of sequences

	# Puts fasta file in memory and parses it
	with open(fastafile,"r") as f:
		for seq in SeqIO.parse(f,"fasta"):
			records.append(seq)

	# retrieves upstream sequence of each gene
	for u in uplist:
		seqid = u[-1]
		for rec in records:
			if rec.id == seqid:
				seq = rec.seq[u[0]-1:u[1]] # extract sequence, beware of indexes, as BioPython indexes from 0
				u.append(seq)

	return uplist

def writeFasta(filename,upseqs):
	"""Write a fastafile from upstream sequences list returned by retrievePos."""

	records = []

	for u in upseqs:
		start = u[0]
		end = u[1]
		strand = u[2]
		gene = u[3] # gene name
		seqid = u[4] # name of scaffold from which the gene is extracted
		seq = u[5] # Seq object

		feat = SeqFeature(FeatureLocation(start,end),strand=strand)
		
		if strand == 1:
			strand = "+"
		elif strand == -1:
			strand = "-"
		else:
			strand = "?"

		ident = seqid+"|"+gene+"|"+str(start)+"-"+str(end)+"|"+strand

		rec = SeqRecord(seq,id=ident,name=gene,features=[feat],description="")

		records.append(rec)

	print "Writing file..."
	SeqIO.write(records,filename,"fasta")
	print "File {} written !".format(filename)
