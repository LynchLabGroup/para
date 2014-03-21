#!/bin/bash
# Script for automatic multialignment using MUSCLE

for FILE in `ls para/data/families/WGD1/CDS/`;
do
FILEBASENAME=$(echo $FILE | cut -d. -f1)
./muscle3.8.31_i86darwin64 -in para/data/families/WGD1/CDS/$FILE -out para/results/multialign/$FILEBASENAMEalign.fasta
done
