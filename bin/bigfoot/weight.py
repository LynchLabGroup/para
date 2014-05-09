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
	def __init__(self,line):
		"""
		line - Pass a sequence line of a StatAlign file
		"""
		self._name = '' #name of the sequence
		self._seq = '' #aligned sequence
		self._w = [] #list of [nt, phyloscore, alignscore] index show position

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

	def weight(self,pred,mpd):
		"""
		Create a list of each base of the sequence object with its prediction and alignment scores.

		pred -- file with phylogenetic prediction scores
		mpd -- file with all alignment scores (.mpd format)
		"""
		weighted = []

		# Create a list of each base with its prediction weight
		with open(pred,"r") as f:
			for c in self._seq:
				num = 0.0
				if c != "-":
					num = f.readline().rstrip("\n") #advances line
					num = float(num) #convert string line to float coefficient
				weighted.append([c,num])

		# Add alignment score
		with open(mpd,"r") as f:
			j = 0
			for i,line in enumerate(f.readlines()):
				try:
					float(line)
				except ValueError:
					j += 1 #count number of lines without number
			f.seek(0) #rewinds the file
			for i,line in enumerate(f.readlines()):
				if i>=j:
					weighted[i-j].insert(2,float(line)) #insert the alignment score

		self._w = weighted

	def motifs(self, thre, size, align):
		"""
		If sequence has been weighted, shows all positions where prediction is higher than threshold.
		Returns a list of motifs with their starting and ending position in the aligned sequence of given minimum size.

		thre -- phylogenetic prediction score threshold, motifs with score over that would be selected
		size -- minimum motifs size to be returned
		align -- alignment score threshold: sequence need to have a score over threshold to be selected
		"""
		if self._w == []:
			return "You first have to weight the sequence using prediction scores !"
		else:
			pos = 0 #position marker
			realpos = 1 # real position makrer in sequence (not including "-" chars)
			known = [] #list of lists position where weight is higher than threshold [start, end, avg score]
			
			#this loop forces to go through the entire sequence
			while pos<len(self._w):
				
				curr_base = self._w[pos][0]  # base
				curr_phyl = self._w[pos][1]  # phylogenetic score
				curr_ali = self._w[pos][2]  # alignment score

				succ_bases = 0  # counter of successive bases
				tot_phyl = 0.0  # variable to calculate average phylogenetic score
				tot_ali = 0.0  # variable for average alignment score
				wide = [] # simple

				# Loop to identify motifs
				while curr_phyl>thre and curr_ali>align:  # if the position have weight and align scores over thresholds
					if succ_bases == 0:
						wide.append(pos+1)  # add the starting position (real position in the sequence)
						realstart = realpos
					succ_bases += 1
					pos += 1
					curr_base = self._w[pos][0] #base

					if curr_base != "-":
						realpos +=1

					curr_phyl = self._w[pos][1]  # phylogenetic score
					curr_ali = self._w[pos][2]  # alignment score
					tot_phyl += curr_phyl
					tot_ali += curr_ali

				#once the loop identified a motif
				if succ_bases < size: #if it is smaller than expected
					pass
				elif succ_bases > 0: #else if we identified a motif
					wide.append(pos) #adding the exact end position of the loop
					avg = tot_phyl/succ_bases
					alavg = tot_ali/succ_bases
					wide.append(avg) #add average score
					wide.append(alavg) #add alignment score
					wide.append(realstart) # add real position start
					known.append(wide)
				if curr_base != "-":
					realpos +=1
				pos += 1

			#extract the motifs sequences using start and end position and insert them in the returned list
			for extract in known:
				if extract[0] == extract[1]:
					sub = self[extract[0]].upper()
				else:
					sub = self[extract[0]:extract[1]+1].upper() #extract the sequence of interest
				extract.insert(0,sub) #insert in the first position the extracted sequence
			return known

	def get_real_pos(self, pos):
		"""
		Returns the realposition of sequence deleting "-" chars with pos a list index
		"""
		real_seq = self._seq[:pos]
		return len(real_seq.translate(None,"-"))



				



