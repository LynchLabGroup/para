#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that from a gene id of GFF file and specified length, gives you the number to extract the upstream sequence of this gene

### IMPORTS ###
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

def retrieve_up(geneids,gff_file,fasta_file,length=100):
	"""
	DESCRIPTION:

	Take a list of gene ids, and retrieve the upstream sequence (default = 100nt) from fasta_file, verifying in gff_file that there is no overlap. Return a list of genes, locations and sequences.

	USAGE:

	geneids -- list of gene ids in gff_file
	gff_file -- a GFF file
	fasta_file -- a Fasta file
	length (optional) -- length of the upstream sequence to extract
	"""

	# GFF file parser, puts GFF file in memory
	rec = load_gff(gff_file) # list of records of each entry in gff file

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
		start = i[0] # pythonic beginning of the gene (begins with position 0)
		end = i[1] # pythonic ending
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
			seqlen = fasta_len(fasta_file,seqid)
			if extract >= fasta_len : #if extract is over the sequence
				extract = seqlen
			for pos in postable:
				if pos[-1] == seqid and pos[0] in range(end+1,extract+1):
					extract = pos[0]
			upstream.append([end,extract,strand,gene,seqid])

	print "Retrieving sequences..."
	upstream = retrieve_seq(fasta_file,upstream)
	
	return upstream

def fasta_len(fasta_file,seqid):
	"""Return length of sequence seqid in fasta_file."""

	length = 0

	# Find the length in fasta_file of seqid
	record = load_fasta(fasta_file)		
	for r in record:
		if record.id == seqid:
			length = len(record.seq) #length of the sequence

	if length == 0:
		print "Sequence {} was not found in {}".format(seqid,fasta_file)
	return length


def retrieve_seq(fasta_file,uplist):
	"""Take retrieve_up return and return the same list with sequences appended."""

	# Puts fasta file in memory and parses it
	records = load_fasta(fasta_file) # list of sequences


	# retrieves upstream sequence of each gene
	for u in uplist:
		seqid = u[-1]
		for rec in records:
			if rec.id == seqid:
				seq = rec.seq[u[0]:u[1]] # extract sequence, beware of indexes, as BioPython indexes from 0
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

def extract_cds(fasta_file,gff_file,output=None):
	"""Returns a list of Coding Sequences extracted from gff_file and fasta_file."""
	if output ==None:
		output = "cds.fa"
	# Load gff file in memory
	rec = load_gff(gff_file)

	fasta = load_fasta(fasta_file)

	# Create a list
	cds = retrieve_pos("CDS",rec) # list of CDS

	# Append sequences at the end of each entry in cds
	for c in cds:
		start = c[0]
		end = c[1]
		seq_id = c[-1]
		seq = seq_extract(seq_id,start,end,fasta)
		c.append(seq)

	
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
		while cds[i+1][4] == gene_id and i+1 < len(cds)-1:
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

	write_fasta(output,genes)

	return cds,genes

def seq_extract(seq_name,start,end,fasta_rec):
	"""Retrieves the DNA sequence in seq_name with positions start and end. fasta_rec is the results of load_fasta"""
	seq = Seq("")
	for f in fasta_rec:
		if f.id == seq_name:
			seq = f.seq[start:end]
	return seq

def retrieve_pos(seq_type,gff_rec):
	"""Return a list of positions sequences of given type using a parsed gff."""

	positions = []
	for r in gff_rec:
		seq_name = r.id # Scaffold name
		feat = r.features # list of all genes in scaffold
		for f in feat:
			gene_name = f.id
			subfeat = f.sub_features # list of all mRNA from this gene
			for s in subfeat:
				subsubfeat = s.sub_features # list of all exons, CDSs and introns
				for ss in subsubfeat:
					if ss.type == seq_type:
						start = ss.location.start.position # Beware it uses position in sequence using pythonic indexes
						end = ss.location.end.position
						strand = ss.location.strand
						phase = int("".join(ss.qualifiers["phase"]))

						posinfo = [start,end,strand,phase,gene_name,seq_name]
						positions.append(posinfo)
	return positions


def load_gff(gff_file):
	"""Returns a list of parsed gff."""
	with open(gff_file,"r") as f:
		print "Parsing file {}...".format(gff_file)
		rec = []
		for line in GFF.parse(f):
			rec.append(line)

	return rec

def load_fasta(fasta_file):
	"""Returns a list of parsed fasta."""
	with open(fasta_file,"r") as f:
		seqs = []
		for seq in SeqIO.parse(f,"fasta"):
			seqs.append(seq)
	return seqs