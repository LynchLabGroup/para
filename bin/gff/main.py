#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gff

def main():

	#family_file = raw_input("Family file? ")
	
	family_file = "data/families/tet_bi_sex_caud_orthoparalogons_WGD1.txt"
	
	# Load family_file in memory
	fam_parser = gff.family_parse(family_file,True)

	length = len(fam_parser[0]) - 1
	print "Detected {} columns apart from name".format(length)

	files = []

	# fasta_files = raw_input("Enter separated by spaces, the names of each fasta file separated by spaces according to header\n")
	# fasta_files = fasta_files.split(" ")

	# gff_files = raw_input("Enter separated by spaces, the names of each gff file separated by spaces according to header\n")
	# gff_files = gff_files.split(" ")

	fasta_files = ["data/tetraurelia/ptetraurelia_mac_51.fa","data/tetraurelia/ptetraurelia_mac_51.fa","data/biaurelia/biaurelia_V1-4_assembly_v1.fasta","data/biaurelia/biaurelia_V1-4_assembly_v1.fasta","data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta","data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta","data/caudatum/caudatum_43c3d_assembly_v1.fasta"]

	gff_files = ["data/tetraurelia/tetraurelia51_EuGene_annotation.gff3","data/tetraurelia/tetraurelia51_EuGene_annotation.gff3","data/biaurelia/biaurelia_V1-4_annotation_v1.gff3","data/biaurelia/biaurelia_V1-4_annotation_v1.gff3","data/sexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3","data/sexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3","data/caudatum/caudatum_43c3d_annotation_v1.gff3"]
	
	# load all fasta files in memory
	fasta_rec = []
	for f in fasta_files:
		fasta_rec.append(gff.load_fasta(f)) 

	# load all gff files in memory
	gff_rec = []
	for g in gff_files:
		gff_rec.append(gff.load_gff(g)) 

	# retrieve all CDS from all gff files and load them in memory
	cds_rec = []
	for g in gff_rec:
		cds_rec.append(gff.retrieve_pos("CDS",g))
	
	# extract all genes for each family
	for fam in fam_parser:
		genes = []
		for i,gene_name in enumerate(fam):
			if i < len(fam)-1 and gene_name != ".":
				genes.append(gff.extract_cds(fasta_rec[i],gff_rec[i],gene_name,cds_rec[i]))
			elif i == len(fam)-1:
				name = fam[i]

		gff.write_fasta("data/"+name+".fasta",genes)


	# genes = extract_cds(fasta_rec,gff_rec,None,cds)

	# write_fasta(output,genes)

if __name__ == "__main__":
	main()
