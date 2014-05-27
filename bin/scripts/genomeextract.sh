#!/bin/bash
# Extract all found motifs in all species
# USAGE:
# genomeextract.sh

extract(){
    local tabfile=$1
    local name=$2
    local gff=$3
    local assembly=$4

    /N/u/jgout/Quarry/bin/extractGenesWithExactMotif.fast $gff GFF Parent $assembly 250 $tabfile > $name.res 2> $name.log

}

LOCATION="results/motifs_extraction/"
FILE="results/26may14uniquememe.txt"
NAME="26may14MEMEMotifs"

extract $FILE $LOCATION$NAME.bi data/biaurelia/biaurelia_cds.gff data/biaurelia/biaurelia_V1-4_assembly_v1.fasta
#sed -e 's/./G/6' $LOCATION$NAME.bi.res > $LOCATION$NAME.bi.res

extract $FILE $LOCATION$NAME.ca data/caudatum/caudatum_cds.gff data/caudatum/caudatum_43c3d_assembly_v1.fasta
#sed -e 's/./G/6' $LOCATION$NAME.ca.res > $LOCATION$NAME.ca.res

extract $FILE $LOCATION$NAME.sex data/sexaurelia/sexaurelia_cds.gff data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta
#sed -e 's/./G/7' $LOCATION$NAME.sex.res > $LOCATION$NAME.sex.res

extract $FILE $LOCATION$NAME.tet data/tetraurelia/tetraurelia_cds.gff data/tetraurelia/ptetraurelia_mac_51.fa
#sed -e 's/./G/8' $LOCATION$NAME.tet.res > $LOCATION$NAME.tet.res

# Code for BigFoot motifs
# <(awk '{print $1$2, $9}' results/15may14consensuscomp.txt)
