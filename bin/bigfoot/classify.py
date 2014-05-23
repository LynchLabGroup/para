#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Contains functions to compute GC content and length of motifs using a list
of sequences.
USAGE:
python classify.py input_file output_file [-f]
"""

from __future__ import print_function

import argparse
from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACDNA
from Bio.SeqUtils import GC123
from Bio.Seq import Seq
from itertools import product


def gc_length(input_file, output, field=None, header=None):
    """
    Given an input tab-delimited field, retrieve sequences at the given field
    (default = last field of tab-delimited file) computes GC content and length
    and outputs file. Header specifies if input has header
    """
    gapped = Gapped(ExtendedIUPACDNA(), '-')  # Extend alphabet to allow gaps

    if field is None:
        field = -1
    if header is None:
        header = True

    with open(input_file, "r") as in_file:
        sequences = []
        number = 0
        for line in in_file:
            if not header:
                line = line.split("\t")
                sequences.append(Seq(line[field], gapped))
            if header and number > 0:
                line = line.split("\t")
                sequences.append(Seq(line[field], gapped))
            number += 1

    with open(output, "w") as o_file:

        dinuc_list = ["".join(x) for x in product("ATCG", repeat=2)]
        header = "AoverAT\tGCm\tlength"
        for nuc in "ATCG":
            header += "\t{}.freq".format(nuc)
        for dinuc in dinuc_list:
            header += "\t{}.freq".format(dinuc)
        header += "\tSeq"

        print(header, file=o_file)
        for seq in sequences:
            gc = GC123(seq)
            at = AoverAT(seq)
            nuc_freq = nucleotide(seq)
            dinuc_freq = dinucleotide(seq)

            line = "{}\t{}\t{}".format(at, gc[0], len(seq)-1)
            for nuc in "ATCG":
                line += "\t{}".format(float(nuc_freq[nuc]))
            for dinuc in dinuc_list:
                line += "\t{}".format(float(dinuc_freq[dinuc]))
            line += "\t{}".format(seq.tostring())

            line = line.rstrip("\n")
            print(line, file=o_file)


def AoverAT(sequence):
    """
    Return the A over AT ratio in sequence
    """
    AT = 0
    A = 0
    for base in sequence:
        if base == "A":
            A += 1
            AT += 1
        elif base == "T":
            AT += 1

    return float(A)/float(AT)


def nucleotide(sequence):
    """Return a dictionary of nucleotide frequency"""
    nuc_list = ["A", "T", "G", "C"]
    nuc_freq = {}
    for nuc in nuc_list:
        count = sequence.count(nuc)
        nuc_freq[nuc] = float(count)/len(sequence)

    return nuc_freq


def dinucleotide(sequence):
    """Returns a dictionary of frequencies of dinucleotide in the sequence."""

    dinuc_list = ["".join(x) for x in product("ATCG", repeat=2)]
    dinuc_freq = {}

    for dinuc in dinuc_list:
        count = sequence.count(dinuc)
        dinuc_freq[dinuc] = float(count)/(len(sequence)/2)

    return dinuc_freq


def main():
    """Main program, containing the command-line parser"""

    # ## Command-line Parser ###
    parser = argparse.ArgumentParser(description="Command-line program to\
        compare motifs outputs of MEME and BigFoot.")

    parser.add_argument("input_file", help="Input tab-delimited file with\
        sequences")
    parser.add_argument("output", help="Output file with seq, GC and length")
    parser.add_argument("-f", "--field", type=int, help="Field number of \
        sequences", default=None)
    parser.add_argument("-head", "--header", help="Specifies if file had header",
                        action="store_true")

    args = parser.parse_args()

    gc_length(args.input_file, args.output, args.field, args.header)

if __name__ == "__main__":
    main()