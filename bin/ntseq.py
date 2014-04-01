#!/usr/bin/env python
# -*- coding: utf-8 -*-
import gff.parse_gff_v2 as pg
import gff.family as fam

def main():
	"""Program to extract all CDSs."""
	
	family_file = "data/families/tet_bi_sex_caud_orthoparalogons_WGD2.txt"
	
	header = True
	num = 7
	translate = False
	
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
	fam.family_cds(fam_parser,gff_rec,fasta_rec,"data/families/WGD2/CDS/nt/",translate=translate)


if __name__ == "__main__":
	main()