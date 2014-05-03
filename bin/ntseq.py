#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import gff.parse_gff_v2 as pg
import gff.family as fam

def main():
	"""Program to extract all CDSs."""
	parser = argparse.ArgumentParser(description="Setup program to extract Coding Sequences of given family")

	parser.add_argument("-f", "--fam", help="specify the minimum number of family members needed. (default: %(default)s)", type=int, default=7)
	parser.add_argument("-t", "--trans", help="If included retrieved sequences are translated.", action="store_true")
	parser.add_argument("-head", "--header", help="specify if family file has a header", action="store_true")
	parser.add_argument("-loc", "--location", help="the folder where to retrieve CDSs (default: %(default)s)", default="data/families/WGD2/CDS/nt/")

	family_file = "data/families/tet_bi_sex_caud_orthoparalogons_WGD2.txt"

	args = parser.parse_args()
	
	header = args.header
	num = args.fam
	translate = args.trans
	loc = args.location
	
	print "Parameters used: {} genes at least.".format(num)	

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
		print "New file..."
		fasta_rec[fk] = pg.load_fasta(fasta_files[fk])
	print "Done."

	# load all gff files in memory
	gff_rec = {}
	print "Loading GFF files..."
	for gk in gff_files.keys():
		print "New file..."
		gff_rec[gk] = pg.load_gff(gff_files[gk],["CDS"]) 
	print "Done."

	# extract all CDSs
	fam.family_cds(fam_parser,gff_rec,fasta_rec,loc,translate=translate)


if __name__ == "__main__":
	main()