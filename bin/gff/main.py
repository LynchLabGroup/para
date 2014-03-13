#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gff
import family as fam

def main():

	#family_file = raw_input("Family file? ")
	
	family_file = "data/families/tet_bi_sex_caud_orthoparalogons_WGD2.txt"
	
	# Load family_file in memory
	args = {}
	args["header"] = True
	args["num"] = 7
	fam_parser,spec = fam.family_parse(family_file,**args)

	print "{} species detected".format(spec)
	fasta_files = {}

	fasta_files["PSEX"] = "data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta"
	fasta_files["PTET"] = "data/tetraurelia/ptetraurelia_mac_51.fa"
	fasta_files["PCAU"] = "data/caudatum/caudatum_43c3d_assembly_v1.fasta"
	fasta_files["PBI"] = "data/biaurelia/biaurelia_V1-4_assembly_v1.fasta"

	# fasta_files = raw_input("Enter separated by spaces, the names of each fasta file separated by spaces according to header\n")
	# fasta_files = fasta_files.split(" ")

	# gff_files = raw_input("Enter separated by spaces, the names of each gff file separated by spaces according to header\n")
	# gff_files = gff_files.split(" ")

	gff_files = {}
	
	gff_files["PSEX"] = "data/sexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3"
	gff_files["PTET"] = "data/tetraurelia/tetraurelia51_EuGene_annotation.gff3"
	gff_files["PCAU"] = "data/caudatum/caudatum_43c3d_annotation_v1.gff3"
	gff_files["PBI"] = "data/biaurelia/biaurelia_V1-4_annotation_v1.gff3"
	
	# load all fasta files in memory
	fasta_rec = {}
	for fk in fasta_files.keys():
		fasta_rec[fk] = gff.load_fasta(fasta_files[fk])

	# load all gff files in memory
	gff_rec = {}
	for gk in gff_files.keys():
		gff_rec[gk] = gff.load_gff(gff_files[gk]) 

	# retrieve all CDS from all gff files and load them in memory
	cds_rec = {}
	for gk in gff_rec.keys():
		cds_rec[gk]=gff.retrieve_pos("CDS",gff_rec[gk])
	
	# extract all genes for each family
	# fam.family_cds(fam_parser,fasta_rec,gff_rec,cds_rec,"data/families/WGD1")

	# extract upstream sequences
	fam.family_upstream(fam_parser,fasta_rec,gff_rec,400,"data/families/WGD2/upstream")


	# genes = extract_cds(fasta_rec,gff_rec,None,cds)

	# write_fasta(output,genes)

if __name__ == "__main__":
	main()
