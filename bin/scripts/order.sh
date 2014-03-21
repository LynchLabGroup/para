#!/bin/bash
# Script for conserving order

for FILE in `ls results/multialign/align*`;
do
# Extract central name of file without prefix nor extension
FILEBASENAME=$(echo $FILE | cut -d. -f1 | cut -d- -f2)
FASTA=$(echo $FILE | cut -d- -f2)
python bin/stable.py data/families/WGD1/CDS/$FASTA $FILE > results/multialign/stable-$FILEBASENAME.fasta
done