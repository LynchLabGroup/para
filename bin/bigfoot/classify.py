#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Contains functions to compute GC content and length of motifs using a list
of sequences.
USAGE:
python classify.py input_file output_file [-f]
"""

import argparse
from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import ExtendedIUPACDNA
from Bio.SeqUtils import GC
from Bio.Seq import Seq


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
        o_file.write("Seq\tGC\tlength")
        for seq in sequences:
            line = "{}\t{}\t{}\n".format(seq.tostring(), GC(seq), len(seq))


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
    parser.add_argument("-h", "--header", help="Specifies if file had header",
                        action="store_true")

    args = parser.parse_args()

    gc_length(args.input_file, args.output, args.field, args.header)

if __name__ == "__main__":
    main()