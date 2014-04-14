#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This scripts makes a dictionnary from interesting files

def dict_block(filename):
	"""
	Return a dictionary of blocks with list of numbers of motifs sizes
	"""

	blocks = {}
	with open(filename, "r") as f:
		for line in f.readlines():
			if line.startswith("WGD"):
				name = line.rstrip("\n")
				blocks[name] = []
			else:
				s = line.rstrip("\n")
				blocks[name].append(int(s))
	return blocks

def int_fam(filename,size_thre):
	"""
	From a file returns a dictionary of motifs sizes per family using size size_thre
	"""
	g = dict_block(filename)

	copy = dict(g)

	for k, sizes in g.items():
		copy[k] = []
		for item in sizes:
			if item >= size_thre:
				copy[k].append(item)

		if copy[k] == []:
			del copy[k]

	return copy