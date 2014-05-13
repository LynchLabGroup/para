#!/bin/bash
# Script to erase folders of discarded families
# USAGE: (be sure to be in "para" folder!)
# erasediscard.sh 

awk '{print "results/"$1}' data/families/WGD2/upstream/discard.txt | xargs rm -r
