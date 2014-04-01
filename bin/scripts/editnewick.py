#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re

def main(input_tree,output_tree):

	# Read tree and create a string out of it
	with open(input_tree,"r") as f:
		l = []
		for line in f.readlines():
			l.append(line)

	l = "".join(l)

	# Edit tree
	newtree = re.sub("\)[0-9].[0-9]*.\:","):",l)

	# Write output tree

	with open(output_tree,"w") as g:
		g.write(newtree)
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Simple program to edit newick tree to be usable with BigFoot.\n Attention tree has to be in Newick format, see http://evolution.genetics.washington.edu/phylip/newicktree.html for specifications.")

	parser.add_argument("i", help="input newick tree")
	parser.add_argument("o",help="output newick tree")

	args = parser.parse_args()

	main(args.i,args.out)
