#!/bin/bash
# Tidy up folder after translatorX script
# USAGE
# tidyup targetlocation

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d

TARGET=$1

find $TARGET -name "*.nt[1-3]*_ali.fasta" -delete
find $TARGET -name "*.aa*" -delete
find $TARGET -name "*.html" -delete
find $TARGET -name "*.muscle.log" -delete