#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# FAMILY FILE -> all functions and files related to gene familes
#
### IMPORTS ###
import gff

### CLASSES ###

class GeneFamily():
	"""
	Object with gene names, and name of the family
	"""

	def __init__(self,line):
		"""Takes a line of a gene family file: a tab-delimited object with all the genes from all species, and first occurrence of WGD is transformed into family name"""
		
		self._name = ""
		self._genes = []

		line = line.rstrip("\n").rstrip("\r").split("\t")

		name_index = first_substring(line,"WGD") # Assumes that first occurrence of "WGD"

		self._name = line[name_index]

		for s in line:
			if "WGD" not in s and s != ".":
				self._genes.append(s)

	def __str__(self):
		return "Family Name: {}\nGenes: {}".format(self._name,self._genes)

	def __repr__(self):
		return self.__str__()

	def __getitem__(self,index):
		return self._genes[index]

	def __getslice__(self,start,end):
		return [self._genes[index] for index in range(start,end)]

	def __len__(self):
		return len(self._genes)

	def __iter__(self):
		for g in self.genes():
			yield g

	def name(self):
		"""Return name of the family"""
		return self._name

	def genes(self):
		"""Return the list of genes"""
		return self._genes

	def species(self):
		"""Return set of species code in family."""

		spec = set()

		for n in self.genes():
			spec.add(get_species(n))

		return spec

	def name_change(self,new_name):
		"""Allow to change family name."""
		print "Former name: {}".format(self._name)
		self._name = new_name
		print "New name: {}".format(self.name()) 

### FUNCTIONS ###

def first_substring(the_list,substring):
	"""Return index of first occurrence of substring in list else return -1."""
	for i,s in enumerate(the_list):
		if substring in s:
			return i
	return -1

def get_species(string):
	"""Get species code of a string. Returns -1 if string length is less than 4"""

	l = len(string)
	if l > 4:
		string = string.split("G")
		species = string[0]
		
		letters = [c for c in species if c.isalpha()] # take only letters

		species = "".join(letters)

		if len(species) > 4:
			species = species[0:4] # Returns only the first four characters
		return species
	else:
		return -1

def family_parse(family_file,num=None,header=None):
	"""Takes a tab-delimited file which gene families on each line.
	By default, assume the file has a header. Last column of the file should be gene family name.
	Return a list of each column where family has at least num number of genes. Return list of species found in file
	"""
	if num == None:
		num = 0

	if header == None:
		header = True

	header_line = []
	with open(family_file,"r") as f:
		family_list = []
		for i,line in enumerate(f.readlines()):
			if header == True and i == 0:
				header_line = GeneFamily(line)
				print "Detected species: {}".format(header_line.species())
				spec = header_line.species()
			
			if header == False or i != 0:
				family = GeneFamily(line)
				if len(family) >= num: # If the family have at least num genes
					family_list.append(family)

	return family_list,spec

def family_cds(family_list,fasta_rec_list,gff_rec_list,cds_rec_list,location):
	"""Take a list of families, with associated dict of fasta_rec, gff_rec and cds_rec. And write Fasta sequences of all CDSs of each family in precised location."""
	for fam in family_list:
		genes = []
		for i,gene in enumerate(fam):
			spec = get_species(gene)
			gene_extract = gff.extract_cds(fasta_rec_list[spec],gff_rec_[spec],gene,cds_rec_list[spec]) # Return a list of list
			genes.append(gene_extract[0]) # extract a simple list
		
		gff.write_fasta(location+fam.name()+".fasta",genes)

def family_upstream(family_list,fasta_rec_list,gff_rec_list,length,location):
	"""Take a list of families, with associated list of fasta_rec and gff_rec. And write upstream Fasta sequences of each family in precised location using length."""

	for fam in family_list:
		up = []
		for i,gene in enumerate(fam):
			spec = get_species(gene)
			gene_list = [gene]
			up_extract = gff.retrieve_up(gene_list,gff_rec_list[spec],fasta_rec_list[spec],length) # Return a list of list
			up.append(up_extract[0]) # extract a simple list

		gff.write_fasta(location+fam.name()+".fasta",up)


