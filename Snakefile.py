# -*- coding: utf-8 -*-
# Snakefile.py file for Snakemake
import os
import glob

RES = [os.path.splitext(famfile)[0]+".comp" for famfile in glob.glob("results/WGD2ANC0000*/WGD2ANC0[0-9][0-9][0-9][0-9].motifs")]
EVAL = 0.001
PREDTHRE = 0.9
ALITHRE = 0.8

rule all:
    input: RES

rule comp:
    """Rule to produce motifs comparison between Bigfoot and MEME motifs"""
    input: "bin/bigfoot/memecomp.py", "{family}.motifs", "{family}.meme.motifs"
    output: "{family}.comp"
    shell: "python2 {input} -e {EVAL} -o {output}"

rule bigfoot_motifs:
    """Parse BigFoot's outputs to produce .motifs files"""
    input: "bin/bigfoot/setup.py", "{family}.fasta.mpd", "{family}.fasta.pred"
    output: "{family}.motifs"
    shell "python2 {input} -o {output} -t {PREDTHRE} -a {ALITHRE}"