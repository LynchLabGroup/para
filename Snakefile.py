# -*- coding: utf-8 -*-
# Snakefile.py file for Snakemake
import os
import glob

RES = [os.path.splitext(famfile)[0]+".comp" for famfile in glob.glob("results/WGD2ANC0000*/WGD2ANC0[0-9][0-9][0-9][0-9].motifs")]
EVAL = 0.001
PREDTHRE = 0.9
ALITHRE = 0.8
MINWIDTH = 4
NMOTIFS = 5
BFPARAM = "10000,20000,1000"

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
    shell: "python2 {input} -o {output} -t {PREDTHRE} -a {ALITHRE}"

rule meme_motifs:
    """Find motifs in sequence with MEME."""
    input: "{family}.fasta"
    output: "{family}.meme.motifs"
    shell: "meme -minw {MINWIDTH} -nmotifs {NMOTIFS} -evt {EVAL} -dna -text {input} >> {output}"

rule bigfoot:
    """Use BigFoot on given files"""
    input: "{family}.fasta", "{family}.newick"
    output: "{family}.fasta.mpd", "{family}.fasta.pred"
    shell: "java -jar ../BigFoot/BigFoot.jar -t {input[1]} -p={BFPARAM} {input[0]}"

rule gene_tree:
    """Create a new tree to be compatible with BigFoot"""
    input: "bin/scripts/editnewick.py", "{family}CDS.phyl_phyml_tree.txt"
    output: "{family}.newick"
    shell: "python2 {input} {output}"

rule phyml_tree:
    """Compute ML gene tree."""
    input: "{family}CDS.phyl"
    output: "{family}CDS.phyl_phyml_tree"
    shell: "phyml -i {input}"

rule fasta_to_phylip:
    """Convert Fasta to Phylip sequences."""
    input: "{family}CDS.e.fasta"
    output: "{family}CDS.phyl"
    shell: "perl bin/scripts/ConvertFastatoPhylip.pl {input} {output}"

rule match_seqs:
    """Matching sequences between CDS and family."""
    input: "bin/scripts/extractmatch.py", "{family}.fasta", "{family}CDS.fasta"
    output: "{family}CDS.e.fasta"
    shell "python2 {input} {output}"

rule transform_CDS_header:
    """Reduce CDS header."""
    input: "bin/scripts/fastaheader.py", "{family}CDS.nt_ali.fasta"
    output: "{family}CDS.fasta"
    shell: "python2 {input} '|' {output}"

rule align_CDS:
    """Aligning CDSs"""
