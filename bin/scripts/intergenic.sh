#!/bin/bash
# Script for computing intergenic distances
python bin/gff/intergenicDistrib.py data/biaurelia/biaurelia_V1-4_annotation_v1.gff3 data/biaurelia/biaurelia_V1-4_assembly_v1.fasta PBI > results/intergenic_pbi.txt

python bin/gff/intergenicDistrib.py data/caudatum/caudatum_43c3d_annotation_v1.gff3 data/caudatum/caudatum_43c3d_assembly_v1.fasta PCAU > results/intergenic_pcau.txt

python bin/gff/intergenicDistrib.py data/sexaurelia/sexaurelia_AZ8-4_annotation_v1.gff3 data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta PSEX > results/intergenic_psex.txt

python bin/gff/intergenicDistrib.py data/tetraurelia/ptetraurelia_CDS_v1.gff3 data/tetraurelia/ptetraurelia_assembly_v1.fa PTET > results/intergenic_ptet.txt

cat results/intergenic* > results/intergenic.csv
