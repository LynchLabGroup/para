# -*- coding: utf-8 -*-
# Snakefile.py file for Snakemake
import os
import glob

RES = [os.path.splitext(famfile)[0]+".comp" for famfile in glob.glob("results/WGD2ANC0000*/WGD2ANC0[0-9][0-9][0-9][0-9].motifs")]
BASENAME = [os.path.basename(family).rstrip(".comp")+".fasta" for family in RES]
UPSTREAM = "data/families/WGD2/upstream"
CDSLOC = "data/families/WGD2/CDS/nt"
RESULTS = "results"
MINLEN = 15
MAXLEN = 250
MINFAMMEMBER = 4
LOWFAMMEMBER = 4
EVAL = 0.001
PREDTHRE = 0.9
ALITHRE = 0.8
MINWIDTH = 4
NMOTIFS = 5
BFPARAM = "10000,20000,1000"

rule all:
    input: RES

rule meme_lengths:
    input: all
    threads: 1
    output: "results/MEMEOutputsLengths.txt"
    shell: "bin/scripts/memelen.sh > {input}"

rule int_motifs:
    input: all
    threads: 1
    output: "{RESULTS}/BFMotifsSummary.txt"
    shell: "cd results/ && echo 'Looking for interesting motifs...'' && \
      .. /bin/scripts/intmotifs.sh BFMotifsSummary.txt && cd .."
rule comp:
    """Rule to produce motifs comparison between Bigfoot and MEME motifs"""
    input: "{family}.motifs", "{family}.meme.motifs"
    output: "{family}.comp"
    shell: "python2 bin/bigfoot/memecomp.py {input} -e {EVAL} -o {output}"

rule bigfoot_motifs:
    """Parse BigFoot's outputs to produce .motifs files"""
    input:  "{family}.fasta.mpd", "{family}.fasta.pred"
    output: "{family}.motifs"
    shell: "python2 bin/bigfoot/setup.py {input} -o {output} -t {PREDTHRE} -a {ALITHRE}"

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
    input: "{family}.CDS.phyl_phyml_tree.txt"
    output: "{family}.CDS.newick"
    shell: "python2 bin/scripts/editnewick.py {input} {output}"

rule phyml_tree:
    """Compute ML gene tree."""
    input: "{family}.CDS.phyl"
    output: "{family}.CDS.phyl_phyml_tree"
    shell: "phyml -i {input}"

rule fasta_to_phylip:
    """Convert Fasta to Phylip sequences."""
    input: "{family}.CDS.e.fasta"
    output: "{family}.CDS.phyl"
    shell: "perl bin/scripts/ConvertFastatoPhylip.pl {input} {output}"

rule match_seqs:
    """Matching sequences between CDS and family."""
    input: "{family}.fasta", "{family}.CDS.fasta"
    output: "{family}.CDS.e.fasta"
    shell: "python2 bin/scripts/extractmatch.py {input} {output}"

rule transform_CDS_header:
    """Reduce CDS header."""
    input: "{family}.CDS.nt_ali.fasta"
    output: "{family}.CDS.fasta"
    shell: "python2 bin/scripts/fastaheader.py {input} '|' {output}"

rule align_CDS:
    """Aligning CDSs"""
    input: "{family}.fasta", "{family}.CDS.fasta"
    params: prefix="{family}.CDS"
    output: "{family}.CDS.nt_ali.fasta"
    shell: "perl bin/scripts/translatorx_vLocal.pl -c 6 -i {input} -o {params.prefix}"

rule make_dir:
    """Make given dir"""
    input: UPSTREAM+"/{family}.fasta"
    output: RESULTS+"/{family}/"
    shell: "mkdir -p {output}"

rule edit_upstream:
    """Edit upstream sequences and move them into results folder."""
    input: UPSTREAM+"/{family}.fasta", make_dir
    output: RESULTS+"/{family}/{family}.fasta"
    shell: "python2 bin/scripts/fastaheader.py {input[0]} '|' {output}"



rule retrieve_CDS:
    """Retrieve CDSs according to families"""
    threads: 1
    shell: "python2 bin/ntseq.py -f {MINFAMMEMBER} --header -loc {CDSLOC}/"

rule retrieve_up:
    """Retrieve upstream sequences."""
    threads: 1
    shell: "python2 bin/gff/main.py -l {MAXLEN} -ml {MINLEN} -f {LOWFAMMEMBER} -mf {MINFAMMEMBER} -loc {UPSTREAM}/ --head"
