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

	def __getitem__(self,num):
		"""
		Returns position num in  object's sequence. Num should be the "natural" position of nucleotide in sequence, beginning with 1.
		"""

		pos = self._seq[num-1]
		return pos

	def __getslice__(self,beg,end):
		"""
		Returns a slice of object's sequence with "natural" numbers.

		'ATTTGC'[2:6] = 'TTTG'
		"""

		extract = self._seq[beg-1:end]
		return extract

	def name(self):
		"""
		Returns name of sequence
		"""
		return self._name

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

	def motifs(self,thre,size):
		"""
		If sequence has been weighted, shows all positions where prediction is higher than threshold.
		Returns a list of motifs with there starting and ending position in the aligned sequence of given minimum size
		"""
		if self._w == []:
			return "You first have to weight the sequence using prediction scores !"
		else:
			pos = 0 #position marker
			known = [] #list of lists position where weight is higher than threshold [start, end, avg score]
			
			#this loop forces to go through the entire sequence
			while pos<len(self._w):
				nt = self._w[pos][0]
				w = self._w[pos][1]
				t = 0
				score = 0.0
				wide = []
				while w>thre: #if the position appears to have a weight higher than threshold
					if t==0:
						wide.append(pos+1) #add the starting position (in term of real position in the sequence)
					t += 1
					pos += 1
					nt = self._w[pos][0]
					w = self._w[pos][1]
					score += w
				
				if t < size:
					pass
				elif t > 0:
					wide.append(pos) #adding the exact end position of the loop
					avg = score/t
					wide.append(avg) #add average score
					known.append(wide)
				pos += 1

			#extract the motifs sequences using start and end position and insert them in the returned list
			for e in known:
				if e[0] == e[1]:
					sub = self[e[0]]
				else:
					sub = self[e[0]:e[1]+1] #extract the sequence of interest
				e.insert(0,sub)
			print "Identified motifs: {}".format(known)
			return known



				



