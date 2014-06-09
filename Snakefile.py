# -*- coding: utf-8 -*-
# Snakefile.py file for Snakemake
import os
import glob
RES = [os.path.splitext(famfile)[0]+".comp" for famfile in glob.glob("results/WGD2ANC0000*/WGD2ANC0[0-9][0-9][0-9][0-9].motifs")]

rule all:
    input: RES

rule comp:
    input: "{family}.motifs", "{family}.meme.motifs"
    output: "{family}.comp"
    shell: "echo {input} && echo {output}"
