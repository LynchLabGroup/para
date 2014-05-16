#!/bin/bash
# Extract all found motifs in all species
# USAGE:
# genomeextract.sh

extract(){
    local motif=$1
    local name=$2
    local gff=$3
    local assembly=$4

    /N/u/jgout/Quarry/bin/extractGenesWithExactMotif $gff GFF Parent $assembly 250 $motif > $name.res 2> $name.log

}

while read -r line; 
do
    fields=(${line//\/})  # Splitting line in array
    echo ${fields[0]}
    # To solve non existing files problem

    extract ${fields[1]} results/motifs_extraction/${fields[0]}bi data/biaurelia/biaurelia_V1-4_annotation_v1.gff3 data/biaurelia/biaurelia_V1-4_assembly_v1.fasta
done < <(awk '{print $1$2, $9}' results/15may14consensuscomp.txt | head)
