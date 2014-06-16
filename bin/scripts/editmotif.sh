#!/bin/bash
# Script to edit motif file
# Usage:
# editmotif.sh sourcefile destination

motiffile=$1
destination=$2
echo "Writing new motif file $destination"
(echo -e "A\tC\tG\tT" && awk '{if(NR>1) print $0}' $motiffile) > $destination
echo "Done."