#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gff

def main():

	fasta = raw_input("What is the path to the fasta file? ")
	gff = raw_input("What is the path to the gff file? ")
	output = raw_input("What should be the output file? ")

	gff_rec = gff.load_gff(gff)
	fasta_rec = gff.load_fasta(fasta)

	cds = gff.retrieve_pos("CDS",gff_rec)

	genes = extract_cds(fasta_rec,gff_rec,None,cds)

	write_fasta(output,genes)

if __name__ == "__main__":
	main()