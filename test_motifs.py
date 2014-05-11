#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Test code to manage modifications for motif search

from bin.bigfoot import parser
s = 'results/WGD2ANC00051/WGD2ANC00051.fasta.'
p = parser.SeqParser(s+".mpd",s+".pred")
p.parse()
p._parse[0].weight(s+".pred",s+".mpd")
