#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 17:47:23 2014

@author: Rekyt
"""


class WeightSeq(object):
    """
    Class for parsed StatAlign alignment sequence lines with three attributes:

    name -- name of the sequence
    seq -- aligned sequence, capitals show determined motidfs
    w -- list of pairs of nt and their phylogenetic footpringtin coefficient
    """
    def __init__(self, line):
        """
        line - Pass a sequence line of a StatAlign file
        """
        self._name = ''  # name of the sequence
        self._seq = ''  # aligned sequence
        self._w = []  # list of [nt, phyloscore, alignscore] index show pos

        # Parse line
        line = line.split("\t")
        name = line[0].strip(" ")  # obtain name of sequence

        seq = line[1].strip("\n")  # aligned sequence, capitals are motifs

        # Initialize object
        self._name = name
        self._seq = seq

    def __str__(self):
        """
        Show name, length and sequence of sequence
        """
        return "name: {}\nlength: {}\nseq: {}".format(self._name,
                                                      len(self._seq), self._seq)

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, num):
        """
        Returns position num in  object's sequence. Num should be the "natural"
        position of nucleotide in sequence, beginning with 1.
        """

        pos = self._seq[num-1]
        return pos

    def __getslice__(self, beg, end):
        """
        Returns a slice of object's sequence with "natural" numbers.

        >>>s = WeightSeq("TESTSEQ\tATTTGC")
        >>>s[2:6]
        'TTTG'
        """

        extract = self._seq[beg-1:end]
        return extract

    def rawseq(self):

        return self._seq.translate(None, "-").upper()

    def name(self):
        """
        Returns name of sequence
        >>>s = WeightSeq("TESTSEQ\tATTTGC")
        >>>s.name
        "TESTSEQ"
        """
        return self._name

    def weight(self, mpd, pred):
        """
        Create a list of each base of the sequence object with its prediction
        and alignment scores.

        mpd -- file with all alignment scores (.mpd format)
        pred -- file with phylogenetic prediction scores
        """
        assert mpd.endswith(".mpd"), pred.endswith(".pred")
        weighted = []

        # Create a list of each base with its prediction weight
        with open(pred, "r") as f:
            for c in self._seq:
                num = 0.0
                if c != "-":
                    num = f.readline().rstrip("\n")  # advances line
                    num = float(num)  # convert string line to float coeff
                weighted.append([c, num])

        # Add alignment score
        with open(mpd, "r") as f:
            j = 0
            for i, line in enumerate(f.readlines()):
                try:
                    float(line)
                except ValueError:
                    j += 1  # count number of lines without number
            f.seek(0)  # rewinds the file
            for i, line in enumerate(f.readlines()):
                if i >= j:
                    weighted[i-j].insert(2, float(line))  # insert ali score

        self._w = weighted

    def motifs(self, thre, align):
        """
        If sequence has been weighted, shows all positions where prediction is
        higher than threshold.
        Returns a list of motifs with their starting and ending position in
        the aligned sequence of given minimum size.

        thre -- phylogenetic prediction score threshold, motifs with score over
                that would be selected
        size -- minimum motifs size to be returned
        align -- alignment score threshold: seq need to have a score
                 over threshold to be selected
        """
        max_window_size = 8
        if self._w == []:
            return "First weight the sequence using prediction scores!"
        else:
            pos = 0  # position marker
            realpos = 1  # real position marker in seq(not including "-" chars)
            known = []  # list of found motifs

            # this loop forces to go through the entire sequence
            while pos < len(self._w):

                wide = []  # keep track of the motif data

                window_size = 1
                good_scores_pos = 0

                # Look through given window_size bases
                while pos < len(self._w) and window_size < max_window_size:
                    if window_size == 1:
                        realstart = realpos

                    curr_base = self._w[pos][0]  # base
                    curr_phyl = self._w[pos][1]  # phylogenetic score
                    curr_ali = self._w[pos][2]  # alignment score
                    window_size += 1
                    if curr_phyl >= thre and curr_ali >= align:
                        good_scores_pos += 1
                    if curr_base != "-":
                        realpos += 1
                    pos += 1

                # If more than 60% of bases are correct in the window
                if float(good_scores_pos)/max_window_size >= 0.75:
                    start_pos = pos - max_window_size + 1  # Start position

                    # Extend motif to the left side (beginning of the sequence)
                    while start_pos > 0 and self._w[start_pos-1][1] >= thre \
                     and self._w[start_pos+1][2] >= align:

                        start_pos -= 1

                    # Extanding the window of the motif to the right
                    while pos < len(self._w) - 1 and self._w[pos+1][1] >= thre and \
                     self._w[pos+1][2] >= align:

                        if self._w[pos+1] != "-":
                            realpos += 1
                        pos += 1

                    avg_phyl = float(sum(y[1] for y in self._w[start_pos:pos+1])) / \
                        (pos - start_pos + 1)
                    avg_ali = float(sum(y[2] for y in self._w[start_pos:pos+1])) / \
                        (pos - start_pos + 1)

                    wide.append(start_pos+1)
                    wide.append(pos)
                    wide.append(avg_phyl)
                    wide.append(avg_ali)
                    wide.append(realstart)
                    known.append(wide)

                    # Go back from the beginning, sliding by one nt
                    pos = realstart+1
                    realpos = pos + 1
                else:
                    pos = pos - window_size + 2
                    realpos = pos + 1

            # Selection for best hit motif without overlap
            if known != []:
                # sort by phylo score
                known.sort(key=lambda motif: motif[2], reverse=True)

                best = best_hits(known)

                # Loop to extract sequences
                for extract in best:
                    if extract[0] == extract[1]:
                        sub = self[extract[0]].upper()
                    else:
                        sub = self[extract[0]:extract[1]+1].upper()
                    extract.insert(0, sub)  # trailing sequence
            return best

    def get_real_pos(self, pos):
        """
        Returns the realposition of sequence deleting "-" chars
        with pos as list index
        """
        real_seq = self._seq[:pos]
        return len(real_seq.translate(None, "-"))


def best_hits(s_list):
    """
    Return the best hits of motifs without overlapping.
    Need to have a sorted list of len > 0

    >>> l =[[191, 200, 0.8181818181818182, 1.0, 192], \
    [191, 200, 0.8181818181818182, 1.0, 194], \
    [212, 218, 0.8125, 1.0, 212], [190, 200, 0.75, 1.0, 190], \
    [213, 220, 0.7222222222222222, 1.0, 214], \
    [191, 202, 0.6923076923076923, 1.0, 196]]
    >>> best_hits(l)
    [[191, 200, 0.8181818181818182, 1.0, 192], [212, 218, 0.8125, 1.0, 212]]
    """
    from memecomp import overlap

    assert len(s_list) > 0, sorted(s_list, key=lambda x: x[2]) == s_list

    best = []
    best.append(s_list[0])

    for i in range(1, len(s_list)):
        coord = (s_list[i][0], s_list[i][1])
        indexes = set()
        for selected in best:
            s_coord = (selected[0], selected[1])
            index = overlap(s_coord, coord, False)
            indexes.add(index)

        if True not in indexes:
            best.append(s_list[i])

    return best