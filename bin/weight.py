#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 17:47:23 2014

@author: Rekyt
"""

class WeightSeq(object):
	"""
	Class for parsed StatAlign alignment sequence lines with three attributes:

	name - name of the sequence
	seq - aligned sequence, capitals show determined motidfs
	w - list of pairs of nt and their phylogenetic footpringtin coefficient
	"""
	def __init__(self,line):
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

		#Initialize object
		self._name = name
		self._seq = seq

	def __str__(self):
		"""
		Show name, length and sequence of sequence
		"""
		return "name: {}\nlength: {}\nseq: {}".format(self._name,len(self._seq),self._seq)

	def __repr__(self):
		return self.__str__()

	def weight(self,pred):
		"""
		Function to weight the sequence according to pred scores in pred file. The pred file must be a column of each pred score.
		"""
		weighted = []
		with open(pred,"r") as f:
			for c in self._seq:
				num = 0.0
				if c != "-":
					num = f.readline().rstrip("\n") #advances line
					num = float(num) #convert string line to float coefficient
				weighted.append([c,num])

		self._w = weighted

	def motifs(self,thre):
		"""
		If sequence has been weighted, shows all positions where prediction is higher than threshold
		"""
		if self._w == []:
			return "You first have to weight the sequence using prediction scores !"
		else:
			for ind,k in enumerate(self._w):
				nt = k[0] #nucleotide
				w = k[1] #phylogenetic prediction
				if w > thre:
					print "index: {} base: {} mpd: {}".format(ind+1,nt,w)

# pos = 0
# known = [] #list of lists position where weight is higher than threshold [start,end]
# while pos<len(self._w):
# 	nt = self._w[pos][0]
# 	w = self._w[pos][1]
# 	t = 0
# 	wide = []
# 	while w>0.1:
# 		if t==0:
# 			wide.append(pos+1)
# 		t += 1
# 		print "{} time {} pos".format(t,pos+1)
# 		pos += 1
# 		nt = self._w[pos][0]
# 		w = self._w[pos][1]
# 	if t > 0:
# 		wide.append(pos)
# 		known.append(wide)
# 	pos += 1

# for e in known:
# 	if e[0] == e[1]:
# 		sub = self._seq[e[0]-1]
# 	else:
# 		sub = self._seq[e[0]-1:e[1]]
# 	e.insert(0,sub)



				



