#!/bin/bash
# Script to use meme on every families 
# USAGE:
# memecomp.sh
for FAM in `ls results/WGD2ANC*/*fasta | grep "WGD2ANC[0-9]\{5\}\.fasta" | cut -d. -f1`;
do
if [ ! -f $FAM.meme.motifs ]
then
	meme $FAM.fasta -minw 4 -nmotifs 5 -dna -text> $FAM.meme.motifs
fi
python bin/bigfoot/memecomp.py $FAM.motifs $FAM.meme.motifs -t 0.9 -e 0.001 -o $FAM.comp
done
