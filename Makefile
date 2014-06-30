# May 12th 2014 pipeline file
# This pipeline compute for all given families
#
#
include config.mak

.PHONY: all makedir retrieve_upstream retrieve_CDS int_motifs meme_length
all:
	echo "Retrieved_FAM"
	echo $(RETRIEVED_FAM)
	echo "SUBDIRS"
	echo $(SUBDIRS)
	echo "FAM"
	echo $(FAM)
	echo "LOCS"
	echo $(LOCS)
	echo "Motifs"
	echo $(MOTIFS)
# #.NOTPARALLEL: retrieve_upstream retrieve_CDS

# all: $(MOTIFS) $(addsuffix .comp, $(LOCS)) int_motifs meme_length

# $(addsuffix .comp, $(LOCS)): %.comp: bin/bigfoot/memecomp.py %.motifs %.meme.motifs
# 	@echo "Comparing MEME and BF"
# 	@python $^ -e $(EVALUE) -o $@
# 	@echo "Done"

# meme_length:
# 	@echo "Comparing MEME outputs lengths"
# 	@bin/scripts/memelen.sh > results/MEMEOutputsLengths.txt
# 	@echo "Done"

# # Underline Interesting motifs
# int_motifs:
# 	cd $(RESULTS) && echo "Looking for interesting motifs..." && \
# 	  ../bin/scripts/intmotifs.sh BFMotifsSummary.txt &&
# 	 cd ..

# # Parse BigFoot's output
# $(MOTIFS): %.motifs : bin/bigfoot/setup.py %.fasta.mpd
# 	@echo "Parsing bigfoot's output"
# 	@python $^ $(subst mpd,pred, $(word 2, $^)) -o $@ -t $(PRED_THRE) -a $(ALI_THRE)
# 	@echo "Done."

# # Compute Motifs using MEME
# $(addsuffix .meme.motifs, $(LOCS)): %.meme.motifs : %.fasta
# 	@echo "Using MEME.."
# 	@meme -minw $(MEME_MIN_W) -nmotifs $(MEME_N_MOTIFS) -dna -text $^ >> $@
# 	@echo "Done."

# # Compute Motifs using BigFoot
# $(MPD) : %.fasta.mpd : %.fasta %.newick
# 	@echo "Computing Motifs using BigFoot"
# 	@java -jar $(BIGFOOT) -t $(word 2, $^) -p=$(BF_PARAMS) $<
# 	@echo "Done."

# # Operations to obtain .newick tree
# $(addsuffix .newick, $(LOCS)):%.newick: bin/scripts/editnewick.py %.CDS.phyl_phyml_tree.txt
# 	@echo "Editing phylo tree"
# 	@python $^ $@
# 	@echo "Done."
# # Obtaining Phyml tree
# $(addsuffix .CDS.phyl_phyml_tree.txt, $(LOCS)): %.CDS.phyl_phyml_tree.txt: %.CDS.phyl
# 	@echo "Computing ML tree"
# 	@echo $<
# 	@phyml -i $<
# 	@echo "Done."

# # Converting fasta to phylip
# $(addsuffix .CDS.phyl, $(LOCS)): %.CDS.phyl: %.CDS.e.fasta
# 	@echo "Converting fasta file to Phylip..."
# 	@perl bin/scripts/ConvertFastatoPhylip.pl $< $@
# 	@echo "Done."

# # Matching sequences between files
# $(addsuffix .CDS.e.fasta, $(LOCS)): %.CDS.e.fasta :bin/scripts/extractmatch.py %.fasta %.CDS.fasta
# 	@echo "Matching sequences"
# 	@python $^ $@
# 	@echo "Done."

# # Transforming CDS header
# $(addsuffix .CDS.fasta, $(LOCS)): %.CDS.fasta : bin/scripts/fastaheader.py %.CDS.nt_ali.fasta
# 	@echo "Transforming CDS header."
# 	@python $^ "|" $@
# 	@echo "Done."
# # Aligning CDSs
# $(addsuffix .CDS.nt_ali.fasta, $(LOCS)): %.CDS.nt_ali.fasta : %.fasta
# 	@echo "Aligning CDSs."
# 	perl $(TRANSLATORX) -c 6 -i $^ -o $(subst .nt_ali.fasta, , $@)
# 	@echo "Done."
# # Transforming fasta header
# $(addsuffix .fasta, $(LOCS)): %.fasta: bin/scripts/fastaheader.py
# 	$(eval UPSTREAM = $(shell echo $@ | sed 's|$(RESULTS)|$(UP)|'))
# 	@echo "Transforming upstream header"
# 	@python $< $(UPSTREAM) "|" $@
# 	@echo "Done."
# # Make appropriate directory
# makedir:
# 	@echo "Making correct dir"
# 	@mkdir -p $(SUBDIRS)

# retrieve_CDS: bin/ntseq.py
# 	@echo "Retrieving CDSs"
# 	@python $^ -f $(FAM_SIZE) --header -loc $(CDS)
# 	@echo "Done."

# retrieve_upstream: bin/gff/main.py
# 	@echo "Retrieving upstream sequences..."
# 	@python $^ -l $(MAX_LEN_UP) -ml $(MIN_LEN_UP) -f $(FAM_SIZE) -mf $(MIN_FAM_SIZE) -loc $(UP) --head
# 	@echo "Done."
