#!/bin/bash
# Script to look at how many interesting motifs have been found
# USAGE:
# intmotifs file_to_write

t=$(for FAM in `ls WGD2ANC0*/*motifs | xargs wc -l | awk '{sub(/.WGD2ANC[0-9]*.motifs/,"")}1 && !/total/' | awk '{if($1 > 2) print $2}'`;
do
	echo $FAM
	awk '{if($11 != "" && $11 != 0.0) print $11}' $FAM/$FAM.motifs
done)

echo "$t" | head
