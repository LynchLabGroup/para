#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test bigfoot module
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'bigfoot'))

import weight
import seqparser

import unittest
import tempfile


class TestClassWeight(unittest.TestCase):
    def setUp(self):
        self.seq = weight.WeightSeq("test\tATTTGC")
        self.gapped = weight.WeightSeq("gapped\tgcg-gggccc----atgcggg")

    def load_gapped(self):
        mpd = "1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n\
            1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n0.0\n0.0\n0.0\n0.0\n0.0\n0.0"
        pred = "0.0\n0.95\n1.0\n1.0\n0.8\n0.95\n1.0\n1.0\n\
            0.95\n0.9\n0.0\n0.0\n0.9\n0.9\n1.0\n1.0\n"
        self.file_mpd = tempfile.NamedTemporaryFile(delete=False)
        self.file_mpd.write(mpd)
        self.file_mpd.close()
        self.file_pred = tempfile.NamedTemporaryFile(delete=False)
        self.file_pred.write(pred)
        self.file_pred.close()

    def test_individual_positions(self):
        self.assertEqual(self.seq[1], "A")
        self.assertEqual(self.seq[5], "G")
        self.assertEqual(self.seq[6], "C")

        with self.assertRaises(IndexError):
            self.seq[7]

    def test_get_slices(self):
        self.assertEqual(self.seq[1:3], "ATT")
        self.assertEqual(self.seq[1:], "ATTTGC")
        self.assertEqual(self.seq[3:1], "")

    def test_rawseq(self):
        self.assertEqual(self.seq.rawseq(), "ATTTGC")
        self.assertEqual(self.gapped.rawseq(), "GCGGGGCCCATGCGGG")

    def test_name(self):
        self.assertEqual(self.seq.name(), "test")

    def test_get_real_pos(self):
        self.assertEqual(self.seq.get_real_pos(1), 1)
        self.assertEqual(self.gapped.get_real_pos(4), 3)

    def test_weight_pred(self):
        self.load_gapped()
        self.gapped.weight(self.file_mpd.name, self.file_pred.name)
        real_pred = [0.0, 0.95, 1.0, 0.0, 1.0, 0.8, 0.95, 1.0, 1.0,
                     0.95, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.9, 0.9,
                     1.0, 1.0]
        self.assertEqual([base[1] for base in self.gapped._w], real_pred)

    def test_weight_mpd(self):
        self.load_gapped()
        self.gapped.weight(self.file_mpd.name, self.file_pred.name)
        real_mpd = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0]
        self.assertEqual([base[2] for base in self.gapped._w], real_mpd)

    def test_motifs_various_thresholds(self):
        self.load_gapped()
        self.gapped.weight(self.file_mpd.name, self.file_pred.name)

        low_threshold = self.gapped.motifs(0.0, 0.0)
        no_motifs = self.gapped.motifs(1.0, 1.0)
        motifs = self.gapped.motifs(0.8, 0.9)

        self.assertListEqual(low_threshold, [['GCG-GGGCCC----ATGCGGG',
            0, 21, 0.5880952380952381, 0.7142857142857143, 14]])
        self.assertListEqual(no_motifs, [])
        self.assertEqual(motifs, [['GGGCCC', 4, 10, 0.9500000000000001, 1.0, 4]])


class TestWeightFunc(unittest.TestCase):
    def test_best_hits_normal(self):
        l = [[3, 8, 0.9], [3, 8, 0.0], [4, 8, 0.7]]
        best = weight.best_hits(l)
        self.assertEqual(best, [[3, 8, 0.9]])

    def test_best_hits_no_list(self):
        self.assertRaises(TypeError, weight.best_hits, [])

    def test_best_hits_score_not_float(self):
        self.assertRaises(TypeError, weight.best_hits, [[0, 0, 2]])

    def test_best_hits_score_too_high(self):
        self.assertRaises(ValueError, weight.best_hits, [[0, 0, 14.0]])

    def test_best_hits_score_too_low(self):
        self.assertRaises(ValueError, weight.best_hits, [[0, 0, -1.0]])


class TestSeqParser(unittest.TestCase):
    """Test class SeqParser from bigfoot/parser.py."""
    def setUp(self):
        mpd = "first\tgcg-gggccc----atgcggg\nsecond\tgcg-gggccc----atgcggg\n\
            1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n\
            1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n1.0\n0.0\n0.0\n0.0\n0.0\n0.0\n0.0"
        pred = "0.0\n0.95\n1.0\n1.0\n0.8\n0.95\n1.0\n1.0\n\
            0.95\n0.9\n0.0\n0.0\n0.9\n0.9\n1.0\n1.0\n"
        self.file_mpd = tempfile.NamedTemporaryFile(delete=False)
        self.file_mpd.write(mpd)
        self.file_mpd.close()
        self.file_pred = tempfile.NamedTemporaryFile(delete=False)
        self.file_pred.write(pred)
        self.file_pred.close()
        print dir(seqparser)
        self.par = seqparser.SeqParser(self.file_mpd, self.file_pred)

    def test_parsing_not_empty(self):
        self.assertNotEqual(self.par.parse(), [])

    def test_parsing_good(self):
        self.assertEqual(self.par.parse(), [])

