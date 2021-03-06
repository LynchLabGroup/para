#!/bin/bash
# Extract all found motifs in all species
# USAGE:
# genomeextract.sh location tabfile prefix
# Location: location where to extract the genes
# Tabfile: simple file with all the motifs to search, one per line
# Prefix: prefix name of the file for each species

extract(){
    local tabfile=$1
    local name=$2
    local gff=$3
    local assembly=$4

    /N/u/jgout/Quarry/bin/extractGenesWithExactMotif.fast $gff GFF Parent $assembly 250 $tabfile > $name.res 2> $name.log

}

replace(){
    local file=$1
    local charnum=$2
    local low=$(($charnum-1))
    local high=$(($charnum+1))


    echo LOW: $low
    echo HIGH: $high
    # Create temporary file with name replaced
    awk '{$2=substr($2,1, '$low') "G" substr($2, '$high'); print}' $file > temp

    mv temp $file
}

LOCATION=$1
FILE=$2
NAME=$3

extract $FILE $LOCATION$NAME.bi data/biaurelia/biaurelia_cds.gff data/biaurelia/biaurelia_V1-4_assembly_v1.fasta
replace $LOCATION$NAME.bi.res 6

extract $FILE $LOCATION$NAME.ca data/caudatum/caudatum_cds.gff data/caudatum/caudatum_43c3d_assembly_v1.fasta
replace $LOCATION$NAME.ca.res 6

extract $FILE $LOCATION$NAME.sex data/sexaurelia/sexaurelia_cds.gff data/sexaurelia/sexaurelia_AZ8-4_assembly_v1.fasta
replace $LOCATION$NAME.sex.res 7

extract $FILE $LOCATION$NAME.tet data/tetraurelia/tetraurelia_cds.gff data/tetraurelia/ptetraurelia_mac_51.fa
replace $LOCATION$NAME.tet.res 8
