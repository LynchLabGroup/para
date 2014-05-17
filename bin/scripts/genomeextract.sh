#!/bin/bash
# Extract all found motifs in all species
# USAGE:
# genomeextract.sh

extract(){
    local motif=$1
    local name=$2
    local gff=$3
    local assembly=$4

    /N/u/jgout/Quarry/bin/extractGenesWithExactMotif $gff GFF Parent $assembly $motif 250 > $name.res 2> $name.log

}

while read -r line; 
do
    fields=(${line//\/})  # Splitting line in array
    echo ${fields[0]}
    # To solve non existing files problem

    extract ${fields[1]} results/motifs_extraction/${fields[0]}bi data/biaurelia/biaurelia_cds.gff data/biaurelia/biaurelia_V1-4_assembly_v1.fasta
    extract ${fields[1]} results/motifs_extraction/${fields[0]}ca data/caudatum/caudatum_cds.gff data/caudatum/caudatum_43c3d_assembly_v1.fasta
    extract ${fields[1]} results/motifs_extraction/${fields[0]}sex data/sexaurelia/sexaurelia_cds.gff data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta
    extract ${fields[1]} results/motifs_extraction/${fields[0]}bi data/tetraurelia/tetraurelia_cds.gff data/tetraurelia/ptetraurelia_mac_51.fa
done < <(awk '{print $1$2, $9}' results/15may14consensuscomp.txt | head)
