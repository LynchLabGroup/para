# May 2nd 2014 pipeline file
# This pipeline compute for all given families
#
#

SUBDIRS = $(shell find results/ -type d -name 'WGD2ANC00002')

FAM = $(subst results, ,$(SUBDIRS))
LOCS = $(join $(SUBDIRS), $(FAM))
MOTIFS = $(addsuffix .motifs, $(LOCS))
MPD = $(addsuffix .fasta.mpd, $(LOCS))
PRED = $(addsuffix .fasta.pred, $(LOCS))


.PHONY: all 

all: $(MOTIFS)

$(MOTIFS): bin/bigfoot/setup.py $(MPD) $(PRED)
	@echo "Parsing bigfoot's output"
	python $^ -o $@ -s 4 -t 0.9 -a 0.8
	@echo "Done."

