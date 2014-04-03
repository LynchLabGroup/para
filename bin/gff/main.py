#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import parse_gff_v2 as pg
import family as fam

def main(length=None, minlength=None, num=None, header=None,location=None):
	"""
	Retrieve all upstream sequences according to a family file, gff files and sequences. Extract upstream sequences with at least minlength nt and max length nt, using family with at least num members. Header tells if the file has a header or not.
	"""
	
	family_file = "data/families/tet_bi_sex_caud_orthoparalogons_WGD2.txt"
	if header == None:
		header = True
	if length == None:
		length = 250
	if minlength == None:
		minlength = 25
	if num == None:
		num = 7
	if location == None:
		location = "data/families/WGD2/upstream/"
	
	print "Parameters used: {} genes at least, retrieve {}nt max, {}nt min".format(num,length, minlength)	

	# Load family_file in memory
	fam_parser,spec = fam.family_parse(family_file,num,header)

	fasta_files = {}

	fasta_files["PSEX"] = "data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta"
	fasta_files["PTET"] = "data/tetraurelia/ptetraurelia_mac_51.fa"
	fasta_files["PCAU"] = "data/caudatum/caudatum_43c3d_assembly_v1.fasta"
	fasta_files["PBI"] = "data/biaurelia/biaurelia_V1-4_assembly_v1.fasta"

	gff_files = {}
	
	gff_files["PSEX"] = "data/sexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3"
	gff_files["PTET"] = "data/tetraurelia/tetraurelia51_EuGene_annotation.gff3"
	gff_files["PCAU"] = "data/caudatum/caudatum_43c3d_annotation_v1.gff3"
	gff_files["PBI"] = "data/biaurelia/biaurelia_V1-4_annotation_v1.gff3"
	
	# load all fasta files in memory
	fasta_rec = {}
	print "Loading fasta files..."
	for fk in fasta_files.keys():
		fasta_rec[fk] = pg.load_fasta(fasta_files[fk])
	print "Done."

	# load all gff files in memory
	gff_rec = {}
	print "Loading GFF files..."
	for gk in gff_files.keys():
		print "New file..."
		gff_rec[gk] = pg.load_gff(gff_files[gk],["CDS","gene"]) 
	print "Done."

	# extract upstream sequences
	fam.family_upstream(fam_parser,gff_rec,fasta_rec,length,minlength,location)



if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="Setup program to extract upstream sequences")

	parser.add_argument("-l", "--length", help="maximum length of extracted upstream sequences. (default: %(default)s)", type=int, default=250)
	parser.add_argument("-ml", "--minlen", help="minimum length of extracted upstream sequences. (default: %(default)s)", type=int, default=1)
	parser.add_argument("-f", "--fam", help="specify the minimum number of family members needed. (default: %(default)s)", type=int, default=7)
	parser.add_argument("--head", help="specify that the file has a header.", action="store_true")
	parser.add_argument("-loc", "--location",help="location of upstream sequences files (default: %(default)s)", default="data/families/WGD2/upstream/")

	args = parser.parse_args()

	main(args.length,args.minlen,args.fam,args.head,args.location)
