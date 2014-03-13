#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=1,vmem=16gb,walltime=20:00
#PBS -M mgrenie@indiana.edu
#PBS -m abe
#PBS -N GeneFamily2
#PBS -j oe

python ~/bin/gff/main.py