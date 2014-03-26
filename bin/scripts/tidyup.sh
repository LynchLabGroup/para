#!/bin/bash
# Tidy up folder after translatorX script
# USAGE
# translatorx 

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d

TARGET="results/translatorx/"

# Look for all files matching pattern in $FILES, and everything that does match the pattern is erased
for FILE in `ls $TARGET | xargs -n1 basename | grep -v ".nt_ali.fasta$"`;
do

echo $FILE | grep "^WGD2ANC*" | xargs rm

done