#!/bin/bash
# Tidy up folder after translatorX script
# USAGE
# translatorx 

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d

TARGET="results/translatorx"

PATTERN="WGD2ANC*.nt_ali.fasta"

echo "Tidying up..."

find $TARGET -not \( -name $PATTERN \) | xargs rm

echo "Done."