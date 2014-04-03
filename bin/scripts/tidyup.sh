#!/bin/bash
# Tidy up folder after translatorX script
# USAGE
# tidyup targetlocation

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d

TARGET=$1

find $TARGET -d 1 -name "*.nt[1-3]*_ali.fasta" -exec rm {} \;
find $TARGET -d 1 -name "*.aa*" -exec rm {} \;
find $TARGET -d 1 -name "*.html" -exec rm {} \;
find $TARGET -d 1 -name "*.muscle.log" -exec rm {} \;