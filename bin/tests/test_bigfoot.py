#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test bigfoot module
import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'bigfoot'))

import weight


class TestWeight(unittest.TestCase):

    def setUp(self):
        self.seq = weight.WeightSeq("test\tATTTGC")
        self.gapped = weight.WeightSeq("gapped\tGCG-GGGCCC")
        # self.real = weight.WeightSeq("")

    def test_individual_positions(self):
        self.assertEqual(self.seq[1], "A")
        self.assertEqual(self.seq[5], "G")
        self.assertEqual(self.seq[6], "C")

        with self.assertRaises(IndexError):
            self.seq[7]

    def test_slices(self):
        self.assertEqual(self.seq[1:3], "ATT")
        self.assertEqual(self.seq[1:], "ATTTGC")
        self.assertEqual(self.seq[3:1], "")

    def test_rawseq(self):
        self.assertEqual(self.seq.rawseq(), "ATTTGC")
        self.assertEqual(self.gapped.rawseq(), "GCGGGGCCC")

    def test_name(self):
        self.assertEqual(self.seq.name(), "test")

    def get_real_pos(self):
        self.assertEqual(self.seq.get_real_pos(1), 1)
        self.assertEqual(self.gapped.get_real_pos(4), 3)

