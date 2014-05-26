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

LOCATION="results/motifs_extraction/MEME/"

while read -r line; 
do
	fields=(${line//\/})  # Splitting line in array
	echo ${fields[0]}
	
	extract ${fields[1]} $LOCATION${fields[0]}bi data/biaurelia/biaurelia_cds.gff data/biaurelia/biaurelia_V1-4_assembly_v1.fasta
	sed -e 's/./G/6' $LOCATION${fields[0]}bi.res > $LOCATION${fields[0]}bi.res
	extract ${fields[1]} $LOCATION${fields[0]}ca data/caudatum/caudatum_cds.gff data/caudatum/caudatum_43c3d_assembly_v1.fasta
	sed -e 's/./G/6' $LOCATION${fields[0]}ca.res > $LOCATION${fields[0]}ca.res
	extract ${fields[1]} $LOCATION${fields[0]}sex data/sexaurelia/sexaurelia_cds.gff data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta
	sed -e 's/./G/7' $LOCATION${fields[0]}sex.res > $LOCATION${fields[0]}sex.res
	extract ${fields[1]} $LOCATION${fields[0]}tet data/tetraurelia/tetraurelia_cds.gff data/tetraurelia/ptetraurelia_mac_51.fa
	sed -e 's/./G/8' $LOCATION${fields[0]}tet.res > $LOCATION${fields[0]}tet.res
done < <(awk '{sub(/\.\/WGD2ANC0[0-9][0-9][0-9][0-9]\/\.\//, "", $6); print}' results/23may14mememotifs.txt | awk '{sub(/\.fasta/, "", $6); print}' | awk '{ if(NR>157) {print $6$5, $4}}')

# Code for BigFoot motifs
# <(awk '{print $1$2, $9}' results/15may14consensuscomp.txt)
