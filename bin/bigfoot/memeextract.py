#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Memeextract.py is a python command-line program to extract relevant TFBS motifs
from MEME XM outputs.
"""


def parse_file(filename, evalue):
    """
    Parse the list of MEME XML outputs and seek motifs with given
    evalue and width.
    """
    sign_motifs = []
    meme_root = load_MEME_xml(filename)
    meme_motifs = meme_root.findall(".//motif")

    # Loop through all motifs in file and extract interesting motifs
    for meme_motif in meme_motifs:
        if float(meme_motif.attrib["e_value"]) <= evalue:
            matrices = find_score_mat(meme_motif)
            for matrix in matrices:
                scores = mat_scores(matrix)
                found_motifs = motifs(scores)

                # Add motif number to list
                found_motifs = [list(motif) for motif in found_motifs]
                [motif.append(meme_motif.attrib["id"])
                 for motif in found_motifs]

                sign_motifs.extend(found_motifs)

    # Add name of data file to list
    sign_motifs = [list(elm) for elm in sign_motifs]
    [elm.append(meme_root.getchildren()[0].attrib["datafile"])
     for elm in sign_motifs]

    return sign_motifs


def load_MEME_xml(filename):
    """Return the root of the parsed XML file."""
    import xml.etree.ElementTree as ET
    xml_tree = ET.parse(filename)
    return xml_tree.getroot()


def find_score_mat(meme_motif):
    """
    From the root of parsed MEME, return list of score elements children of
    probabilities. Assumes that all motifs in file have good scores.
    """
    return meme_motif.findall(".//probabilities/alphabet_matrix")


def mat_scores(alpha_matrix):
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


def best_base(letter_scores):
    """
    From a letter_score dictionary, extract best score. If all letters have
    the same score, returns the first letter in alphabetical order.

    >>> let = {'A':0.4, 'T':0.0, 'C':0.7, 'G':0.9}
    >>> best_base(let)
    'G'
    >>> let = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0}
    >>> best_base(let)
    'A'
    >>> let = {'A':1.0, 'T':1.0, 'C':1.0, 'G':1.0}
    >>> best_base(let)
    'A'
    >>> let = {'A':0.0, 'T':1.0, 'C':1.0, 'G':1.0}
    >>> best_base(let)
    'C'
    """
    letter = max(sorted(letter_scores.keys()), key=(lambda key:
                                                    letter_scores[key]))
    return letter


def best_score(letter_scores):
    """
    Return the score of the best base."""
    letter = best_base(letter_scores)
    return letter_scores[letter]


def motifs(scores, window_size=None, good_pos=None, threshold=None):
    """
    Using a scores list, look for potential motif And return a list of them.
    """
    if window_size is None:
        window_size = 8
    if good_pos is None:
        good_pos = 6
    if threshold is None:
        threshold = 0.8
    pos = 0  # Position in the scores list
    motifs = []
    while pos < len(scores) - window_size + 1:
        valid_pos = 0  # Number of positions with good scores
        curr_mot = ""

        # Sliding window along positions
        for curr_wind in xrange(window_size):
            if curr_wind == 0:
                wind_start = pos  # Record precisely beginning of window
            best_nt = best_base(scores[pos])
            curr_mot += best_nt
            if scores[pos][best_nt] >= threshold:
                valid_pos += 1
            pos += 1

        # What happens if the motif has enough good pos
        if valid_pos >= good_pos:
            start_pos = pos - window_size

            # Extend motif left than right
            while best_score(scores[start_pos-1]) >= threshold \
             and start_pos > 0:
                start_pos -= 1
                curr_mot = best_base(scores[start_pos]) + curr_mot
            while best_score(scores[pos]) >= threshold \
             and pos < len(scores) - 1:
                curr_mot = curr_mot + best_base(scores[pos])
                pos += 1

            # Trim motifs for eventual non relevant bases on the edges
            while best_score(scores[start_pos]) < threshold:
                start_pos += 1
                curr_mot = curr_mot[1:]
            while best_score(scores[pos - 1]) < threshold:
                pos -= 1
                curr_mot = curr_mot[:-1]
            motifs.append((start_pos, pos, curr_mot))
            pos = wind_start + 1
        else:
            pos = pos - window_size + 1

    if motifs != []:
        motifs = best_motifs(motifs, scores)

    return motifs


def best_motifs(motifs_list, scores_list):
    """
    From a list of motifs, return non-overlapping motifs having best scores.
    """
    from memecomp import overlap

    assert type(motifs_list) == list and type(scores_list) == list and\
     motifs_list != [] and scores_list != []

    motifs = [list(motif) for motif in motifs_list]  # Make a mutable copy

    # Compute the score for all motifs
    for motif in motifs:
        motif.insert(2, avg_bestscore(motif[0], motif[1], scores_list))

    # Sort motif by best score
    motifs.sort(key=lambda motif: motif[2], reverse=True)
    best = []
    best.append(tuple(motifs[0]))

    for i in range(1, len(motifs)):
        coord = (motifs[i][0], motifs[i][1])
        indexes = set()
        for selected in best:
            s_coord = (selected[0], selected[1])
            index = overlap(s_coord, coord, False)
            indexes.add(index)

        if True not in indexes:
            best.append(tuple(motifs[i]))
    return best


def avg_bestscore(start, end, scores_list):
    """
    Return the mean of the best scores from start to end included
    """
    assert type(start) == int and type(end) == int and\
    type(scores_list) == list
    score_sum = 0
    for index in xrange(start, end):
        score_sum += best_score(scores_list[index])
    return float(score_sum)/(end - start + 1)


def write_tabfile(motifs_list, out_file):
    """Write a tab-delimited file from a motif list."""
    with open(out_file, "w") as out:
        header = "start\tstop\tavg.score\tmotif\tMEME.motif\tdatafile\n"
        out.write(header)
        for motif in motifs_list:
            line = ""
            for elm in motif:
                line += "{}\t".format(elm)
            line.rstrip("\t")
            line += "\n"

            out.write(line)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Command-line interface to\
        extract motifs from MEME XML outputs.")

    parser.add_argument("evalue", help="motif evalue", type=float)
    parser.add_argument("output_file", help="Output file for motifs")
    parser.add_argument("memefiles", help="One or more meme.xml files",
                        nargs="+")

    args = parser.parse_args()

    int_motifs = []
    # Extend list of found motifs with motif from each files
    for mfile in args.memefiles:
        print mfile
        int_motifs.extend(parse_file(mfile, args.evalue))

    write_tabfile(int_motifs, args.output_file)

if __name__ == "__main__":
    main()







    