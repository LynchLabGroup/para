#!/bin/bash
# Script for diff between different directories
# USAGE:
# diff dir1 dir2

d=$(date "+%Y-%m-%d") #Stores date value into variable
echo $d
echo "Beginning diff..."

out=diff$d.txt

diff -q -r $1 $2 >> $out

echo "Written file $out"

echo "Done."
