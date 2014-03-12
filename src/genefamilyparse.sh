#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=1,walltime=30:00
#PBS -M mgrenie@indiana.edu
#PBS -m abe
#PBS -N GeneFamily
#PBS -j oe

cd ../para/
python ~/bin/gff/main.py