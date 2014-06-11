#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test bigfoot module
import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'bigfoot'))

import weight


class WeightSeq(unittest.TestCase):

    def setup(self):
        self.seq = weight.WeightSeq("test\tATTTGC")
        self.gapped = weight.WeightSeq("gapped\tGCG-GGGCCC")
        self.real = weight.WeightSeq("")

    def test_individual_positions(self):
        self.assertEqual(self.seq[1], "A")
        self.assertEqual(self.seq[5], "G")
        self.assertEqual(self.seq[len(self.seq)], "C")

        with self.assertRaises(IndexError):
            self.seq[len(self.seq)+1]

    def test_slices(self):
        self.assertEqual(self.seq[1:3], "AT")
        self.assertEqual(self.seq[1:], "ATTTGC")
        self.assertEqual(self.seq[3:1], [])

    def test_rawseq(self):
        self.assertEqual(self.seq.rawself.seq(), "ATTTGC")
        self.assertEqual(self.gapped.rawself.seq(), "GCGGGGCCC")

    def test_name(self):
        self.assertEqual(self.seq.name(), "test")

    def get_real_pos(self):
        self.assertEqual(self.seq.get_real_pos(1), 1)
        self.assertEqual(self.gapped.get_real_pos(4), 3)

