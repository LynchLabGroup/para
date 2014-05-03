# May 2nd 2014 pipeline file
# This pipeline compute for all given families
#
#

SUBDIRS = $(addprefix results/, $(shell awk '{print $1}' results/NotRetrieved.txt))

FAM = $(subst results, ,$(SUBDIRS))
LOCS = $(join $(SUBDIRS), $(FAM))
MOTIFS = $(addsuffix .motifs, $(LOCS))
MPD = $(addsuffix .fasta.mpd, $(LOCS))
PRED = $(addsuffix .fasta.pred, $(LOCS))
UP = "data/families/WGD2/upstream/"
UP_S = $(subst upstream/, upstream, $(UP)) 
CDS = "data/families/WGD2/CDS/nt/"


.PHONY: all makedir retrieve_upstream retrieve_CDS

all: retrieve_CDS makedir $(MOTIFS)

# Parse BigFoot's output
$(MOTIFS): bin/bigfoot/setup.py $(MPD) $(PRED)
	@echo "Parsing bigfoot's output"
	python $^ -o $@ -s 4 -t 0.9 -a 0.8
	@echo "Done."

# Compute Motifs using MEME
$(addsuffix .meme.motifs, $(LOCS)): $(addsuffix .fasta, $(LOCS))
	@echo "Using MEME.."
	@meme -minw 4 -nmotifs 5 -dna -text $^ >> $@
	@echo "Done."

# Compute Motifs using BigFoot
$(MPD) $(PRED): $(addsuffix .fasta, $(LOCS)) $(addsuffix .newick, $(LOCS))
	@echo "Computing Motifs using BigFoot"
	@java -jar ../BigFoot/BigFoot.jar -t $(word 2, $^) -p=1000,2000,1000 $<
	@echo "Done."

# Operations to obtain .newick tree
$(addsuffix .newick, $(LOCS)): bin/scripts/editnewick.py $(addsuffix CDS.phyl_phyml_tree.txt, $(LOCS))
	@echo "Editing phylo tree"
	@python $^ $@
	@echo "Done."
# Obtaining Phyml tree
$(addsuffix CDS.phyl_phyml_tree.txt, $(LOCS)): $(addsuffix CDS.phyl, $(LOCS))
	@echo "Computing ML tree"
	@phyml -i $<
	@echo "Done."

# Converting fasta to phylip
$(addsuffix CDS.phyl, $(LOCS)): $(addsuffix CDS.e.fasta, $(LOCS))
	@echo "Converting fasta file to Phylip..."
	@perl bin/scripts/ConvertFastatoPhylip.pl $< $@
	@echo "Done."

# Matching sequences between files
$(addsuffix CDS.e.fasta, $(LOCS)): bin/scripts/extractmatch.py $(addsuffix .fasta, $(LOCS)) $(addsuffix CDS.fasta, $(LOCS))
	@echo "Matching sequences"
	@python $^ $@
	@echo "Done."

# Transforming CDS header
$(addsuffix CDS.fasta, $(LOCS)): bin/scripts/fastaheader.py $(addsuffix CDS.nt_ali.fasta, $(LOCS))
	@echo "Transforming CDS header."
	@python $^ "|" $@
	@echo "Done."
# Aligning CDSs
$(addsuffix CDS.nt_ali.fasta, $(LOCS)): $(addsuffix .fasta, $(LOCS))
	@echo "Aligning CDSs."
	@perl bin/scripts/translatorx_vLocal.pl -c 6 -i $^ -o $(addsuffix CDS, $(LOCS))
	@echo "Done."
# Transforming fasta header
$(addsuffix .fasta, $(LOCS)): bin/scripts/fastaheader.py $(addsuffix .fasta, $(addprefix $(UP_S),$(FAM)))
	@echo "Transforming upstream header"
	@python $^ "|" $@
	@echo "Done."
# Make appropriate directory
makedir:
	@echo "Making correct dir"
	@mkdir -p $(LOCS)

retrieve_upstream: bin/gff/main.py
	@echo "Retrieving upstream sequences..."
	@python $^ -l 250 -ml 15 -f 4 -mf 4 -loc $(UP) --head
	@echo "Done."
retrieve_CDS: bin/ntseq.py
	@echo "Retrieving CDSs"
	@python $^ -f 4 --header -loc $(CDS)
	@echo "Done."
