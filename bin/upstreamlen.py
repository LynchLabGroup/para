#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Program that plot upstream sequences lengths of given files

### IMPORTS ###
from gff import gff_func
from gff import parse_gff_v2 as pg
import csv

### FUNCTIONS ###
def main():
	maxlen = 20000
	output = "results/upstream_lengths.csv"

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

	# Build a dictionnary of all lengths of upstream sequences in all species
	lengths = {}
	for k in gff_rec.keys():
		lengths[k] = gff_func.retrieve_up_len(gff_rec[k],fasta_rec[k],maxlen)
	spec = []
	val = []
	for k in lengths.keys():
		for l in lengths[k]:
			spec.append(k)
			val.append(l)
	
	print "Writing CSV file {}...".format(output)
	with open(output,"w") as f:
		writer = csv.writer(f, delimiter='\t')

		# header row
		writer.writerow(["species","value"])
		for row in len(spec):
			writer.writerow([spec[row],val[row]])
	print "Done."
			
if __name__ == "__main__":
	main()