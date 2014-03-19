#!/usr/bin/env python
# -*- coding: utf-8 -*-
import parse_gff_v2 as pg
import gff_p
import family as fam

def main():

	#family_file = raw_input("Family file? ")
	
	family_file = "data/families/test_families.txt"
	length = 400
	
	
	header = True
	num = 7
	
	print "Parameters used: {} genes at least, retrieve {}nt".format(num,length)	

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
		gff_rec[gk] = pg.load_gff(gff_files[gk]) 
	print "Done."

	# extract upstream sequences
	fam.family_upstream(fam_parser,gff_rec,fasta_rec,400,"data/families/test_families/upstream/")



if __name__ == "__main__":
	main()
