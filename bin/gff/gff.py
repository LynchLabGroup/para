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
	upstream = []
	for g in geneids:
		start = gff_dic[g].start
		end = gff_dic[g].end
		strand = gff_dic[g].strand
		seq_id = gff_dic[g].id
		seq = retrieve_up_single(g,gff_dic,fasta_dic,length)
		upstream.append([start,end,strand,g,seq_id,seq])
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
		
		return fasta_dic[seq_id][extract-1:start-1].reverse_complement().complement() # extract natural position, reverse the return to begin with first base just before the beginning of the gene
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
		start = u[0]
		end = u[1]
		strand = u[2]
		gene = u[3] # gene name
		seqid = u[4] # name of scaffold from which the gene is extracted
		seq = u[-1] # Seq object

		ident = seqid+"|"+gene+"|"+start+"-"+end+"|"+strand

		rec = SeqRecord(seq,id=ident,name=gene,description="")

		records.append(rec)

	print "Writing file..."
	SeqIO.write(records,file_name,"fasta")
	print "File {} written !".format(file_name)

def extract_cds(gff_dic,fasta_dic,gene_name=None):
	"""
	Returns a list of Coding Sequences extracted from gff_dic and fasta_dic.
	
	WARNING: functions strangely with gff files where there are several transcript for a single gene
	"""

	interest = []
	if gene_name == None:
		parents = set()
		for g in gff_dic:
			if gff_dic[g].type == "CDS":
				par = gff_dic[g].attributes["Parent"]
				parents.add(par)

		for p in parents:
			start = gff_dic[p].start
			end = gff_dic[p].end
			strand = gff_dic[p].strand

			# Change name to have gene name
			name = list(p)
			name[-6] = "G"
			name = "".join(name)
			seq_id = gff_dic[p].id
			seq = gp.parse_gff.get_prot_seq(gff_dic,fasta_dic,p)

			interest.append([start,end,strand,name,seq_id,seq])

	else:
		# Extract from a list of names
		print "gene_name = {}".format(gene_name)
		
		for g in gff_dic:
			if gff_dic[g].id == gene_name:
				gene = gff_dic[g]
				start = gene.start
				end = gene.end
				strand = gene.strand
				seq_id = gene.id

				trans = list(gene_name)
				trans[-6] = "T"
				trans = "".join(trans)

				seq = gp.parse_gff.get_prot_seq(gff_dic,fasta_dic,trans)

				interest.append([start,end,strand,name,seq_id,seq])

		if interest == []:
			sys.exit("Gene(s) of interest: {} was not found".format(gene_name))

	return interest

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
		print "No positions were found for sequence of given type."
	else:
		return positions
