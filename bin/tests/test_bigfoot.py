#!/usr/bin/env python
# -*- coding: utf-8 -*-
# test bigfoot module
import unittest
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'bigfoot'))

import weight


class WeightSeq(unittest.TestCase):

    def setUp(self):
        seq = weight.WeightSeq("test\tATTTGC")
    def test_individual_positions(self):
        self.assertEqual(seq[1], "A")
        self.assertEqual(seq[5], "G")
        self.assertEqual(seq[len(seq), "C")

        with self.assertRaises(IndexError):
            seq[len(seq)+1]

    def test_slices(self):
        self.assertEqual(seq[1:3], "AT")
        self.assertEqual(seq[1:], "ATTTGC")
        self.assertEqual(seq[3:1], [])

    def test_rawseq(self):
        self.assertEqual(seq.rawseq(), "ATTTGC")

    def test_name(self):
        self.assertEqual(seq.name(), "test")
        
