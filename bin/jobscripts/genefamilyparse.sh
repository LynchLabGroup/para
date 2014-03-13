#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=2,vmem=14gb,walltime=01:00:00
#PBS -M mgrenie@indiana.edu
#PBS -m abe
#PBS -N GeneFamily
#PBS -j oe

module load python
module load biopython
ls
cd para/
python bin/gff/main.py
