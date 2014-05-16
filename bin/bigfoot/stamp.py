#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to format a sequence alignment file for stamp
"""
from __future__ import print_function
from memecomp import MotifFile
import argparse


def fullpath(fam_name):
    """
    From a family name return the path of the motif file
    """
    local_str = str(fam_name)
    directory = "results/"
    return directory+local_str+"/"+local_str+".motifs"


def extract_ali(input_file, output_file):
    """
    From a allcomp.sh generated file, generate an output file of all alignment
    of all motifs.
    """
    motif_files = {}  # Dictionnary of all loaded motifs file
    with open(input_file, "r") as in_file:
        with open(output_file, "w") as out_file:

            for line in in_file:
                line = line.rstrip("\n").split("\t")

                if line[0] not in motif_files.keys():
                    # Load motif file if not the case
                    motif_files[line[0]] = MotifFile(fullpath(line[0]))
                print(">{}".format(line[0]+line[1]), file=out_file)

                for motif in motif_files[line[0]]:
                    if motif.name == line[1]:
                        break

                # Print the sequence alignment of the motif
                print(motif, file=out_file)


def main():
    """
    Main program function
    """
    # ## Command-line Parser ###
    parser = argparse.ArgumentParser(description="Command-line program to\
        compare generate output for STAMP\
        http://www.benoslab.pitt.edu/stamp/index.php")

    parser.add_argument("allcomp", help="allcomp.sh output file")
    parser.add_argument("out", help="Name of the output file")

    args = parser.parse_args()
    print("Extracting motifs alignments...")
    extract_ali(args.allcomp, args.out)
    print("Done.")

if __name__ == "__main__":
    main()
