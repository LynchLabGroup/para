#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script to load sequences and swap the order of sequences in fasta files
# USAGE:
# python swaporder.py input_fasta output_fasta

from Bio import SeqIO


def load_fasta(input_file):
    """
    Return a list of fasta sequences.
    """
    assert input_file.endswith(".fasta") or input_file.endswith(".fa")
    with open(input_file, "r") as f:
        fasta_list = []
        for rec in SeqIO.parse(f, "fasta"):
            fasta_list.append(rec)

    return fasta_list


def randomize_list(input_list):
    """
    Return a randomized list from input_list.
    """
    import random
    copy = list(input_list)
    random.shuffle(copy)
    return copy


def write_fasta(input_list, output_file):
    """
    Write fasta file from a given list.
    """
    with open(output_file, "w") as output:
        SeqIO.write(input_list, output, "fasta")


def main(input_file, output_file):
    """
    Main program
    """
    print "Loading file {}...".format(input_file)
    in_list = load_fasta(input_file)
    print "Loaded!"

    r_list = randomize_list(in_list)

    print "Writing file {}...".format(output_file)
    write_fasta(r_list, output_file)
    print "Done."

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Simple program to generate\
        random ordered fasta files.")

    parser.add_argument("i", help="input fasta file")
    parser.add_argument("o", help="output fasta file")

    args = parser.parse_args()

    main(args.i,args.o)