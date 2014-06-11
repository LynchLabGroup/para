#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Command-line interface to extract the upstream regions of given list of genes

import gff_func
import parse_gff_v2 as pg
import argparse


def gene_names(inputname):
    with open(inputname, "r") as infile:
        records = [line.rstrip("\n") for line in infile.xreadlines()]
    return records


def extract(inputname, outputname, gff, fasta):
    """Take a list of gene name and extract their upstream regions using
        given gff and fasta files.
    """
    records = gene_names(inputname)
    upstream_seqs = gff_func.retrieve_up(records, gff, fasta, 250, 15)[0]

    # Select only retrieved sequences
    upstream = [seq for seq in upstream_seqs if len(seq) > 1]
    return upstream


def main():
    parser = argparse.ArgumentParser(description="Setup program to extract \
        upstream sequences")

    parser.add_argument("genelist", help="Simple list of genes you want \
        to extract upstream regions, without header")
    parser.add_argument("outputname", help="Fasta output file")
    parser.add_argument("gff", help="GFF file to use")
    parser.add_argument("fasta", help="Genome Fasta file to use")

    args = parser.parse_args()

    gff = pg.load_gff(args.gff, ["CDS", "gene"])
    fasta = pg.load_fasta(args.fasta)
    inputname = args.genelist
    outputname = args.outputname

    upstream = extract(inputname, outputname, gff, fasta)

    gff_func.write_fasta(outputname, upstream)

if __name__ == "__main__":
    main()
