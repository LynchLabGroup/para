# June 30th 2014 pipeline file
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
# Upstream sequences location
UP = data/families/WGD2/upstream/
# CDSs location
CDS = data/families/WGD2/CDS/nt/
# Results directory (where all families will be treated)
RESULTS = results/
# All retrieved family
RETRIEVED_FAM = $(shell find $(UP) -maxdepth 1 -type f -name "WGD2ANC0000*")
# List corresponding subdirectories
SUBDIRS = $(foreach DIR, "${RETRIEVED_FAM[0]}", $(subst $(UP), $(RESULTS), $(basename $(DIR))))
# List all families with just name of families
FAM = $(foreach NAME, $(RETRIEVED_FAM), $(basename $(NAME)))
# Base file names
LOCS = $(join $(SUBDIRS), $(FAM))
# Files
MOTIFS = $(addsuffix .motifs, $(LOCS))
PRED = $(addsuffix .fasta.pred, $(LOCS))
MPD = $(addsuffix .fasta.mpd, $(LOCS))
