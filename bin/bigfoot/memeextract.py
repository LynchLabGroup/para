#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Memeextract.py is a python command-line program to extract relevant TFBS motifs
from MEME XM outputs.
"""


def parse_files(files_list, evalue, minwidth=None, maxwidth=None):
    """
    Parse the list of MEME XML outputs and seek motifs with given
    evalue and width.
    """
    pass


def load_MEME_xml(filename):
    """Return the root of the parsed XML file."""
    import xml.etree.ElementTree as ET
    xml_tree = ET.parse(filename)
    return xml_tree.getroot()


def find_score_matrices(xml_root):
    """
    From the root of parsed MEME, return list of score elements children of
    probabilities.
    """
    assert xml_root.tag == "MEME"
    return xml_root.findall(".//probabilities/alphabet_matrix")


def matrix_scores(alpha_matrix):
    """
    From an alphabet matrix return a data structure containing a list of
    dictonaries containing scores for each letter at each position.
    """
    scores_list = []

    for array in alpha_matrix.iter():
        if array.tag == "alphabet_array":
            letter_scores = {}
            # Loop over each letter in the array ATCG
            for value in array.iter():
                if value.tag == "value":
                    letter = value.attrib["letter_id"][-1]
                    letter_scores[letter] = float(value.text)

            scores_list.append(letter_scores)
    return scores_list

    