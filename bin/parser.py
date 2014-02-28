#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 17:47:23 2014

@author: Rekyt
"""

class WeigthSeq(object):
	"""
	Class for parsed StatAlign alignment sequence lines with three attributes:

	name - name of the sequence
	seq - aligned sequence, capitals show determined motidfs
	w - list of pairs of nt and their phylogenetic footpringtin coefficient
	"""
	def __init__(self,line,pred):
		"""
		line - Pass a sequence line of a StatAlign file
		pred - pass name of pred file
		"""
		self._name = '' #name of the sequence
		self._seq = '' #aligned sequence
		self._w = [] #list of [nt, mpd] index show position

		#Parse line
		line = line.split("\t")
		name = line[0].strip(" ") #obtain name of sequence

		seq = line[1].strip("\n") #aligned sequence, capitals are motifs

		weighted = []
		with open("pred","r") as f:
			for c in seq:
				num = 0.0
				if c != "-":
					num = f.readline().rstrip("\n") #advances line
					num = float(num) #convert string line to float coefficient
				weighted.append([c,num])


		#Initialize object
		self._name = name
		self._seq = seq
		self._w = weighted
