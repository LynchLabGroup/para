# -*- coding: utf-8 -*-
# Snakefile.py file for Snakemake
import os
import glob
RES = [os.path.splitext(famfile)[0]+".comp" for famfile in glob.glob("results/WGD2ANC0000*/WGD2ANC0[0-9][0-9][0-9][0-9].motifs")]
EVAL = 0.001
rule all:
    input: RES

rule comp:
    input: "{family}.motifs", "{family}.meme.motifs"
    output: "{family}.comp"
    shell: "python2 bin/bigfoot/memecomp.py {input} -e {EVAL} -o {output}"
