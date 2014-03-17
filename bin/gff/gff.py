#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that from a gene id of GFF file and specified length, gives you the number to extract the upstream sequence of this gene

### IMPORTS ###
from .. import gff_python as gp # import sibling folder
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

def retrieve_up(geneids,gff_dic,fasta_dic,length=100):
	"""
	DESCRIPTION:

	Take a list of gene ids, and retrieve the upstream sequence (default = 100nt) from fasta_dic, verifying in gff_file that there is no overlap. Return a list of genes, locations and sequences.

	USAGE:

	geneids -- list of gene ids in gff_file
	gff_dic -- a records of GFF file
	fasta_dic -- a Fasta rec
	length (optional) -- length of the upstream sequence to extract
	"""

	# GFF file parser, puts GFF file in memory
	

	interest = [] # list of coordinates of gene of interest
	postable = [] # list of coordinates of all genes
	# Looping into content to identify characteristics of gene of interest
	print "Extracting genes of interest..."
	for k in gff_dic.keys():
		f = gff_dic[k]
		name = f.attributes["ID"]
		start = f.start # retrieve int as position
		end = f.end 
		strand = f.strand

		# Retrieves data from genes of interest
		for g in geneids:
			if name == g:
				interest.append([start,end,strand,g,f.seqid])
		
		# Retrieves all the genes
		if f.type == "gene":
			postable.append([start,end,strand,f.seqid]) #constructing a position matrix

	# Retrieves start and end position to extract, while verifying there is no overlapping genes
	upstream = []
	print "Retrieving upstream locations..."
	for i in interest:
		start = i[0] # pythonic beginning of the gene (begins with position 0)
		end = i[1] # pythonic ending
		strand = i[2] # strand +: plus -: minus
		gene = i[3] # gene id
		seqid = i[4] # sequence id on which gene is located
		extract = 0
		
		if strand == "+":
			extract = start - length # beginning of the extract
			if extract <= 0: # we don't want to extract negative bases!
				extract = 1
			for pos in postable:
				if pos[-1] == seqid and pos[1] in range(extract,start):# if the gene is overlapping on the same scaffold
					extract = pos[1]+1
			upstream.append([extract,start,strand,gene,seqid])# gives the sequence you want to extract upstream of the gene		
		elif strand == "-":
			extract = end + length # end of the extract
			seqlen = fasta_len(fasta_dic,seqid)
			if extract >= seqlen : #if extract is over the sequence
				extract = seqlen
			for pos in postable:
				if pos[-1] == seqid and pos[0] in range(end+1,extract+1):
					extract = pos[0]
			upstream.append([end,extract,strand,gene,seqid])

	print "Retrieving sequences..."
	upstream = retrieve_seq(fasta_dic,upstream)

	# Reverse order of sequences to begin with first base close to gene
	for u in upstream:
		if u[2] == "+":
			u[-1] = u[-1].reverse_complement().complement() #reverse order of sequence
	
	return upstream

def retrieve_up_single(geneid,gff_dic,fasta_dic,length=100):
	"""
	Return upstream sequence of given length, without overlapping for a single gene.
	"""

	for k in gff_dic.keys():
		if k == geneid:
			start = gff_dic[k].start
			end = gff_dic[k].stop
			strand = gff_dic[k].strand
			seq_id = gff_dic[k].id

	if strand == "+":
		extract = start - length # We don't want to include first base of gene
		if extract <= 0:
			extract = 1 #No negative bases !
		for g in gff_dic:
			if g.id == seq_id and g.end in range(extract,start): # g.end is the last position of gene in sequence, whatever its strand
				extract = g.end + 1
			elif g.id == seq_id and g.start in range(extract,start):
				print "Detected overlapping element {}\ntype: {} start: {} end: {}".format(g,g.type,g.start,g.end)
		return fasta_dic[seq_id][extract-1:start-1] # extract natural position
	elif strand == "-":
		extract = end + length
		if extract >= len(fasta_dic[k]):
			extract = len(fasta_dic[k])
		for g in gff_dic:
			if g.id == seq_id and g.start in range(end+1,extract+1):
				extract = g.start
			elif g.id == seq_id and g.end in range(end+1,extract+1):
				print "Detected overlapping element {}\ntype: {} start: {} end: {}".format(g,g.type,g.start,g.end)
		return fasta_dic[seq_id][end:extract]


def fasta_len(fasta_dic,seqid):
	"""Return length of sequence seqid in fasta_rec."""

	length = 0

	# Find the length in fasta_file of seqid	
	for k in fasta_dic.keys():
		if k == seqid:
			length = len(fasta_dic[k]) #length of the sequence

	if length == 0:
		print "Sequence {} was not found".format(seqid)
	return length


def retrieve_seq(fasta_dic,uplist):
	"""Append sequence to the passed list, based on fasta_dic."""

	# Puts fasta file in memory and parses it
	fasta = fasta_dic # list of sequences


	# retrieves upstream sequence of each gene
	for u in uplist:
		seqid = u[-1]
		for k in fasta.keys():
			if k == seqid:
				seq = fasta[k][u[0]-1:u[1]] # extract sequence, beware of indexes, as BioPython indexes from 0
				u.append(seq)

	return uplist

def write_fasta(file_name,upseqs):
	"""Write a fasta_file from upstream sequences list returned by retrieve_up."""

	records = []

	for u in upseqs:
		strand = u[2]
		gene = u[3] # gene name
		seqid = u[4] # name of scaffold from which the gene is extracted
		seq = u[-1] # Seq object

		
		if strand == 1:
			strand = "+"
		elif strand == -1:
			strand = "-"
		else:
			strand = "?"

		ident = seqid+"|"+gene+"|"+strand

		rec = SeqRecord(seq,id=ident,name=gene,description="")

		records.append(rec)

	print "Writing file..."
	SeqIO.write(records,file_name,"fasta")
	print "File {} written !".format(file_name)

def extract_cds(fasta_dic,gff_dic,gene_name=None,cds=None):
	"""Returns a list of Coding Sequences extracted from gff_dic and fasta_dic."""

	fasta = fasta_dic

	# Create a list
	if cds == None:
		cds = retrieve_pos("CDS",gff_dic) # list of CDS

	# Append sequences at the end of each entry in cds
	if gene_name == None:
		for c in cds:
			start = c[0]
			end = c[1]
			seq_id = c[-1]
			seq = seq_extract(seq_id,start,end,fasta)
			c.append(seq)
	else:
		# Extract from a list of names
		print "gene_name = {}".format(gene_name)
		interest = []
		for c in cds:
			if c[4] in gene_name:
				start = c[0]
				end = c[1]
				seq_id = c[-1]
				seq = seq_extract(seq_id,start,end,fasta)
				c.append(seq)
				interest.append(c)
		if interest == []:
			sys.exit("Gene(s) of interest: {} was not found".format(gene_name))
		cds = interest


	
	# This loop assemble translated genes
	genes = []
	i = 0
	while i < len(cds)-1:
		start = cds[i][0]
		strand = cds[i][2]
		phase = cds[i][3]
		gene_id = cds[i][4]
		seq_id = cds[i][5]
		dna_seq = Seq("")
		dna = []
		dna.append(cds[i][-1])


		while i+1 <= len(cds)-1 and cds[i+1][4] == gene_id:
			i += 1
			dna.append(cds[i][-1])

		if strand == 1:
			for d in dna:
				dna_seq += d
		elif strand == -1:
			dna.reverse()
			for d in dna:
				dna_seq += d.reverse_complement()
		
		prot_seq = dna_seq.translate(table=6)

		genes.append([start,phase,strand,gene_id,seq_id,prot_seq])
		i += 1

	return genes

def retrieve_pos(seq_type,gff_dic):
	"""Return a list of positions sequences of given type using a parsed gff."""
	positions = []
	for k in gff_dic.keys():
		if gff_dic[k].type == seq_type:

			## Attributes
			start = gff_dic[k].start
			end = gff_dic[k].end
			strand = gff_dic[k].strand
			phase = gff_dic[k].phase
			name = k
			seq_name = ggf_dic[k].seqid

			## What to append
			posinfo = [start,end,strand,phase,gene_name,seq_name]
			
			positions.append(posinfo)
	if positions == []:
		print "No positions were found. For sequence of given type."
	else:
		return positions
