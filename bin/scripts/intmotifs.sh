#!/bin/bash
# Script to look at how many interesting motifs have been found
# USAGE:
# intmotifs lengthsmotifsfile > intmotifsfile

for FAM in `awk '{if ($1 > 2) print $2}' $1`;
do
	echo $FAM
	awk '{if($11 != "" && $11 != 0.0) print $11}' $FAM/$FAM.motifs
done
