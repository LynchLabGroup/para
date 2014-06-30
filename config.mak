# May 12th 2014 pipeline file
# This config file contains all variables

SUBDIRS = $(shell find results/ -type d -name "WGD2ANC0*") 
FAM = $(subst results, ,$(SUBDIRS))
LOCS = $(join $(SUBDIRS), $(FAM))
MOTIFS = $(addsuffix .motifs, $(LOCS))
MPD = $(addsuffix .fasta.mpd, $(LOCS))
PRED = $(addsuffix .fasta.pred, $(LOCS))
UP = data/families/WGD2/upstream/
CDS = data/families/WGD2/CDS/nt/
EVALUE = 0.001
PRED_THRE = 0.9
ALI_THRE = 0.8
BF_PARAMS = 10000,20000,1000
MEME_MIN_W = 4
MEME_N_MOTIFS = 5
FAM_SIZE = 4
MIN_FAM_SIZE = 4
MIN_LEN_UP = 15
MAX_LEN_UP = 250
