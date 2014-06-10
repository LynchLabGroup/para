#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script to extract highly expressed non
# ribosomal genes -> background

import sys
import os
# Add gff to PYTHONPATH
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'gff'))

import gff_func
import parse_gff_v2 as pg

gff = pg.load_gff("data/tetraurelia/tetraurelia51_EuGene_annotation.gff3", ["CDS","gene"])
fasta = pg.load_fasta("data/tetraurelia/ptetraurelia_mac_51.fa")

with open("results/highnonmatch.txt", "r") as infile:
    records = [line.rstrip("\n") for line in infile.xreadlines()]
upstream_seqs = gff_func.retrieve_up(records, gff, fasta, 250, 15)[0]

# Select only retrieved sequences
upstream = [seq for seq in upstream_seqs if len(seq) > 1]
gff_func.write_fasta("results/ribo.background.fasta", upstream)
