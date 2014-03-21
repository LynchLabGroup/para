#!/bin/bash
# Script for automatic multialignment using MUSCLE
# need to have muscle executable at the root!

for FILE in `ls para/data/families/WGD1/CDS/`;
do
FILEBASENAME=$(echo $FILE | cut -d. -f1)
./muscle3.8.31_i86darwin64 < para/data/families/WGD1/CDS/$FILE > para/results/multialign/align-$FILE -log outalign-$FILEBASENAME
done
