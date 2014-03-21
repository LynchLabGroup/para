#!/bin/bash
# Script for automatic multialignment using MUSCLE

for FILE in `ls para/data/families/WGD1/CDS/`;
do
FILEBASENAME=$(echo $FILE | cut -d. -f1)
./muscle3.8.31_i86darwin64 < para/data/families/WGD1/CDS/$FILE > para/data/families/WGD1/CDS/align-$FILE
done
