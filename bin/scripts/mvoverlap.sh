#!/bin/bash
# Script to move overlapped families

echo "Moving overlapped families..."
for FILE in `cut -f 1 data/families/WGD2/upstream/overlap.txt`;
do

mv data/families/WGD2/upstream/$FILE.fasta data/families/WGD2/upstream/discarded/$FILE.fasta
echo "Moved file $FILE..."
done
echo "Done."