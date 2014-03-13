#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=1,vmem=10gb,walltime=00:20:00
#PBS -M mgrenie@indiana.edu
#PBS -m abe
#PBS -N GeneFamily3
#PBS -j oe

module load python
module load biopython
cd ..
cd para/
python bin/gff/main.py