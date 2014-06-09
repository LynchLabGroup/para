# -*- coding: utf-8 -*-
# Snakemake.py file for Snakemake

RES = "results/"

rule comp:
    input: "{RES}/{family,WGD2ANC0000.*}/{family}.motifs", "{RES}/{family,WGD2ANC0000.*}/{family}.meme.motifs"
    output: "{RES}/{family}/{family}.comp"
    shell: "python bin/bigfoot/memecomp.py {input} -e 0.001 -o {output}"