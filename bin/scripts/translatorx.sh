#!/bin/bash
# Script to automatically use translatorx to align sequences
# USAGE
# translatorx 

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d
echo "Beginning diff..."

SOURCEPATH="data/families/WGD2/CDS/nt/"

TARGET="results/translatorx/"

for FILE in `ls $SOURCEPATHWGD*`
do

FILEBASENAME=$(echo $FILE | cut -d. -f1)

perl bin/translatorx_vLocal.pl -c 6 -i $SOURCEPATH$FILE -o $TARGET$FILEBASENAME
done