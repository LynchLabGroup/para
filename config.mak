# May 12th 2014 pipeline file
# This config file contains all variables

SUBDIRS = $(shell find results/ -type d -name "WGD2ANC00051") 
FAM = $(subst results, ,$(SUBDIRS))
LOCS = $(join $(SUBDIRS), $(FAM))
MOTIFS = $(addsuffix .motifs, $(LOCS))
MPD = $(addsuffix .fasta.mpd, $(LOCS))
PRED = $(addsuffix .fasta.pred, $(LOCS))
UP = data/families/WGD2/upstream/
CDS = data/families/WGD2/CDS/nt/