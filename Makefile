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
UP = data/families/WGD2/upstream/
UP_S = data/families/WGD2/upstream
CDS = data/families/WGD2/CDS/nt/
#OBJ = $(wildcard results/WGD2ANC*/WGD2ANC[0-9]*.motifs)
DIRS = $(wildcard results/WGD2ANC*/)
SRC = $(wildcard $(addsuffix WGD2ANC*.motifs,$(DIRS)))
OBJ = $(filter %.motifs, $(SRC))
.PHONY:  all makedir retrieve_upstream retrieve_CDS

all: $(MOTIFS)
#makedir $(MOTIFS)
#echo $(DIRS)
#	echo $(OBJ)
	# echo $(sources)
	# echo $(objects)
#$(MOTIFS): %.motifs : %.fasta.mpd %.fasta.pred 
#	@echo "Making $@" 

# Parse BigFoot's output
$(MOTIFS): %.motifs : bin/bigfoot/setup.py %.fasta.mpd %.fasta.mpd
	@echo "Parsing bigfoot's output"
	@python $^ -o $@ -s 4 -t 0.9 -a 0.8
	@echo "Done."

# Compute Motifs using MEME
$(addsuffix .meme.motifs, $(LOCS)): %.meme.motifs : %.fasta
	@echo "Using MEME.."
	@meme -minw 4 -nmotifs 5 -dna -text $^ >> $@
	@echo "Done."

# Compute Motifs using BigFoot
$(MPD) :%.fasta.mpd : %.fasta %.newick
	@echo "Computing Motifs using BigFoot"
	@java -jar ../BigFoot/BigFoot.jar -t $(word 2, $^) -p=1000,2000,1000 $<
	@echo "Done."

# Operations to obtain .newick tree
$(addsuffix .newick, $(LOCS)):%.newick: bin/scripts/editnewick.py %CDS.phyl_phyml_tree.txt
	@echo "Editing phylo tree"
	@python $^ $@
	@echo "Done."
# Obtaining Phyml tree
$(addsuffix CDS.phyl_phyml_tree.txt, $(LOCS)): %CDS.phyl_phyml_tree.txt: %CDS.phyl
	@echo "Computing ML tree"
	@echo $<
	@phyml -i $<
	@echo "Done."

# Converting fasta to phylip
$(addsuffix CDS.phyl, $(LOCS)): %CDS.phyl: %CDS.e.fasta
	@echo "Converting fasta file to Phylip..."
	@perl bin/scripts/ConvertFastatoPhylip.pl $< $@
	@echo "Done."

# Matching sequences between files
$(addsuffix CDS.e.fasta, $(LOCS)): %CDS.e.fasta :bin/scripts/extractmatch.py %.fasta %CDS.fasta
	@echo "Matching sequences"
	@python $^ $@
	@echo "Done."

# Transforming CDS header
$(addsuffix CDS.fasta, $(LOCS)): %CDS.fasta : bin/scripts/fastaheader.py %CDS.nt_ali.fasta
	@echo "Transforming CDS header."
	@python $^ "|" $@
	@echo "Done."
# Aligning CDSs
$(addsuffix CDS.nt_ali.fasta, $(LOCS)): %CDS.nt_ali.fasta : %.fasta
	@echo "Aligning CDSs."
	perl bin/scripts/translatorx_vLocal.pl -c 6 -i $^ -o $(subst .nt_ali.fasta, , $@)
	@echo "Done."
# Transforming fasta header
$(addsuffix .fasta, $(LOCS)): %.fasta: bin/scripts/fastaheader.py
	@echo $(UP_S)
	@echo "Transforming upstream header"
	@python $< $(addsuffix .fasta, $(addprefix $(UP_S),$(FAM))) "|" $@
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
