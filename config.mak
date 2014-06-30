# May 12th 2014 pipeline file
# This config file contains all variables

### Parameters ###
# Upstream sequence maximum length
MAX_LEN_UP = 250
# Upstream sequence minimum length
MIN_LEN_UP = 15
# Minimal gene family size to look at the family
FAM_SIZE = 4
# Minimal gene family size (after removing improper genes)
MIN_FAM_SIZE = 4

# Number of motifs to detect with MEME
MEME_N_MOTIFS = 5
# Minimum size of detected motifs with MEME
MEME_MIN_W = 4
# Max evalue of MEME motifs
EVALUE = 0.001

# Parameters for BigFoot (burn-in, cycles, sampling rate)
BF_PARAMS = 10000,20000,1000

# Alignment score threshold to detect motifs
ALI_THRE = 0.8
# Prediction score threshold to detect motifs
PRED_THRE = 0.9

### Files and Directories ###

UP = data/families/WGD2/upstream/
CDS = data/families/WGD2/CDS/nt/
RESULTS = results/
SUBDIRS = $(shell find results/ -type d -name "WGD2ANC0*")
GENES = $(shell find $(UP) -maxdepth 1 -type f -name "WGD2ANC0*")
SUBS = $(subst $(UP), $(RESULTS), $(basename $(GENES)))

MOTIFS = $(addsuffix .motifs, $(LOCS))
LOCS = $(join $(SUBDIRS), $(FAM))
FAM = $(subst results, ,$(SUBDIRS))
PRED = $(addsuffix .fasta.pred, $(LOCS))
MPD = $(addsuffix .fasta.mpd, $(LOCS))
