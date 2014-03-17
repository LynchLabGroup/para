#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ..gff_python import parse_gff as pg
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

	# retrieve all CDS from all gff files and load them in memory
	cds_rec = {}
	for gk in gff_rec.keys():
		cds_rec[gk]=gff.retrieve_pos("CDS",gff_rec[gk])

	# extract upstream sequences
	fam.family_upstream(fam_parser,gff_rec,fasta_rec,400,"data/families/WGD2/upstream/")


if __name__ == "__main__":
	main()
