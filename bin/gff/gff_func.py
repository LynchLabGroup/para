#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that from a gene id of GFF file and specified length, gives you the number to extract the upstream sequence of this gene

### IMPORTS ###
import parse_gff as pg

def retrieve_up(geneids,gff_dic,fasta_dic,length=100):
	"""
	DESCRIPTION:

	Take a list of gene ids, and retrieve the upstream sequence (default = 100nt) from fasta_dic, verifying in gff_file that there is no overlap. Return a list of genes, locations and sequences.

	USAGE:

	geneids -- list of gene ids in gff_file
	gff_dic -- a records of GFF file
	fasta_dic -- a Fasta record
	length (optional) -- length of the upstream sequence to extract
	"""
	upstream = []
	over = []

	# create list of scaffold ids
	scaff = [k for k in gff_dic.keys() if "scaff" in k]
	
	#retrieve upstream sequence of all given genes
	for g in geneids:
		for seq_id in scaff:
			l = find_name(gff_dic,seq_id,g)
			if l != -1:
				start,end,strand = l
				break
		
		seq,overlap = retrieve_up_single(g,seq_id,gff_dic,fasta_dic,length)
		upstream.append([start,end,strand,g,seq_id,seq])
		
		# detect if gene has strange overlap
		if overlap > 0:
			over.append(g)
	return upstream,over 

############################################

def retrieve_up_single(geneid,seqid,gff_dic,fasta_dic,length=100):
	"""
	Return upstream sequence of given length, without overlapping for a single gene.
	"""

	try:
		gff_dic[seqid]
		overlap = 0 #overlapping index
		
		# Find gene coordinates in the scaffold
		start,end,strand = find_name(gff_dic,seqid,geneid)

		# Return gene sequence according to strand as well as number of overlapping genes
		if strand == "+":
			extract = start - length # We don't want to include first base of gene
			if extract <= 0:
				extract = 1 #No negative bases !
			for g in gff_dic[seqid]:
				if g.end in xrange(extract,start): # g.end is the last position of gene in sequence, whatever its strand
					extract = g.end + 1
				elif g.start in xrange(extract,start):
					print "Detected overlapping element {}\ntype: {} start: {} end: {}".format(g,g.type,g.start,g.end)
					overlap += 1
			
			return fasta_dic[seqid][extract-1:start-1].reverse_complement().complement(),overlap # extract natural position, reverse the return to begin with first base just before the beginning of the gene
		elif strand == "-":
			extract = end + length
			if extract >= len(fasta_dic[seqid]):
				extract = len(fasta_dic[seqid])
			for g in gff_dic[seqid]:
				if g.start in xrange(end+1,extract+1):
					extract = g.start
				elif g.end in xrange(end+1,extract+1):
					print "Detected overlapping element {}\ntype: {} start: {} end: {}".format(g,g.type,g.start,g.end)
					overlap += 1
			return fasta_dic[seqid][end:extract],overlap
	except KeyError:
		return -1
	
############################################

def fasta_len(fasta_dic,seqid):
	"""Return length of sequence seqid in fasta_rec."""

	length = 0

	# Find the length in fasta_file of seqid	
	try:
		length = len(fasta_dic[seqid])
		return length
	except KeyError:
		print "Sequence {} was not found".format(seqid)

############################################

def retrieve_seq(fasta_dic,uplist):
	"""Append sequence to the passed list, based on fasta_dic."""

	# Puts fasta file in memory and parses it
	fasta = fasta_dic # list of sequences

	# retrieves upstream sequence of each gene
	for u in uplist:
		seqid = u[-1]
		try:
			seq = fasta[seqid][u[0]-1:u[1]]
			u.append(seq)
		except KeyError:
			print "Key {} was not found.".format(seqid)
	return uplist

############################################

def write_fasta(file_name,upseqs):
	"""Write a fasta_file from upstream sequences list returned by retrieve_up."""
	from Bio.SeqRecord import SeqRecord
	from Bio import SeqIO

	records = []

	for u in upseqs:
		start = u[0]
		end = u[1]
		strand = u[2]
		gene = u[3] # gene name
		seqid = u[4] # name of scaffold from which the gene is extracted
		seq = u[-1] # Seq object

		ident = "{} | {} | {}-{} | {}".format(seqid,gene,str(start),str(end),strand)

		rec = SeqRecord(seq,id=ident,name=gene,description="")

		records.append(rec)

	print "Writing file..."
	SeqIO.write(records,file_name,"fasta")
	print "File {} written !".format(file_name)

############################################

def extract_cds(gff_dic,fasta_dic,par_name=None):
	"""
	Returns a list of Coding Sequences extracted from gff_dic and fasta_dic using parent name.
	WARNING: functions strangely with gff files where there are several transcript for a single gene
	"""

	interest = []
	if par_name == None:
		parents = [k for k in gff_dic.keys() if "scaff" not in k]
		for p in parents:
			
			start = gff_dic[p][0].start
			end = gff_dic[p][-1].end
			strand = gff_dic[p][0].strand


			# Change name to have gene name
			name = list(p)
			name[-6] = "G"
			name = "".join(name)
			seq_id = gff_dic[p][0].seqid
			seq = pg.get_prot_seq(gff_dic,fasta_dic,p)

			interest.append([start,end,strand,name,seq_id,seq])

	else:
		# Extract sequence from a single parent name
		print "par_name = {}".format(par_name)
		
		for k in gff_dic.keys():
			if k == par_name:
				gene = gff_dic[g]
				start = gene.start
				end = gene.end
				strand = gene.strand
				seq_id = gene.seqid
				trans = list(par_name)
				trans[-6] = "T"
				trans = "".join(trans)

				seq = pg.get_prot_seq(gff_dic,fasta_dic,trans)

				interest.append([start,end,strand,name,seq_id,seq])
				break

		if interest == []:
			print "Gene(s) of interest: {} was not found".format(par_name)

	return interest

############################################

def find_name(gff_dic,seqid,geneid):
	"""Return start,end,strand and seqid of gene in seqid."""
	found = False
	for i,r in enumerate(gff_dic[seqid]):
			if r.attributes["ID"] == geneid:
				start = gff_dic[seqid][i].start
				end = gff_dic[seqid][i].end
				strand = gff_dic[seqid][i].strand
				found = True
				break
	if found == False:
		return -1
	else:
		return [start,end,strand]
