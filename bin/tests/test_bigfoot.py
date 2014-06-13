#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test bigfoot module
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'bigfoot'))

import weight
import mock
import unittest
import tempfile
from StringIO import StringIO


class TestWeight(unittest.TestCase):
    def setUp(self):
        self.seq = weight.WeightSeq("test\tATTTGC")
        self.gapped = weight.WeightSeq("gapped\tgcg-gggccc----atgcggg")

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
        self.gapped.weight(self.file_mpd.name, self.file_pred.name)
        real_pred = [0.0, 0.95, 1.0, 0.0, 1.0, 0.8, 0.95, 1.0, 1.0,
                     0.95, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.9, 0.9,
                     1.0, 1.0]
        self.assertEqual([base[1] for base in self.gapped._w], real_pred)

    def test_weight_mpd(self):
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
        self.gapped.weight(self.file_mpd.name, self.file_pred.name)
        real_mpd = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0]
        self.assertEqual([base[2] for base in self.gapped._w], real_mpd)