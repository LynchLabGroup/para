#!/bin/bash
# Script to automatically use translatorx to align sequences
# USAGE
# translatorx 

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d

SOURCEPATH="data/families/WGD2/CDS/nt/"

TARGET="results/translatorx/"

FILES="WGD*"

SOURCE=$SOURCEPATH$FILES


# Look for all files matching pattern in $FILES and use translatorX on them
for FILE in `ls $SOURCE | xargs -n1 basename`;
do

FILEBASENAME=$(echo $FILE | cut -d. -f1)

perl bin/scripts/translatorx_vLocal.pl -c 6 -i $SOURCEPATH$FILE -o $TARGET$FILEBASENAME
done