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

# Parse BigFoot's output
$(MOTIFS): bin/bigfoot/setup.py $(MPD) $(PRED)
	@echo "Parsing bigfoot's output"
	python $^ -o $@ -s 4 -t 0.9 -a 0.8
	@echo "Done."i

# Compute Motifs using MEME
$(addsuffix .meme.motifs, $(LOCS)): $(addsuffix .fasta, $(LOCS))
	meme -minw 4 -nmotifs 5 -dna -text $^ >> $@

# Compute Motifs using BigFoot
$(MPD) $(PRED): $(addsuffix .fasta, $(LOCS)) $(addsuffix .newick, $(LOCS))
	java -jar ../BigFoot/BigFoot.jar -t $(word 2, $^) -p=1000,2000,1000 $<
# Transforming fasta header
$(addsuffix .fasta, $(LOCS)): bin/scripts/fastaheader.py $(addsuffix .fasta, $(addprefix data/families/WGD2/upstream/,$(FAM)))
	python $^ "|" $@

# Operations to obtain .newick tree
$(addsuffix .newick, $(LOCS)): bin/scripts/editnewick.py $(addsuffix CDS.phyl_phyml_tree.txt, $(LOCS))
	python $^ $@

# Obtaining Phyml tree
$(addsuffix CDS.phyl_phyml_tree.txt, $(LOCS)): $(addsuffix CDS.phyl, $(LOCS))
	phyml -i $<

# Converting fasta to phylip
$(addsuffix CDS.phyl, $(LOCS)): $(addsuffix CDS.e.fasta, $(LOCS))
	perl bin/scripts/ConvertFastatoPhylip.pl $< $@

# Matching sequences between files
$(addsuffix CDS.e.fasta, $(LOCS)): bin/scripts/extractmatch.py $(addsuffix .fasta, $(LOCS)) $(addsuffix CDS.fasta, $(LOCS))
	python $^ $@

# Transforming CDS header
$(addsuffix CDS.fasta, $(LOCS)): bin/scripts/fastaheader.py $(addsuffix CDS.nt_ali.fasta, $(LOCS))
	python $^ "|" $@

# Aligning CDSs

$(addsuffix CDS.nt_ali.fasta, $(LOCS)): $(addsuffix .fasta, $(LOCS))
	perl bin/scripts/translatorx_vLocal.pl -c 6 -i $^ -o $(addsuffix CDS, $(LOCS))

retrieve_upstream: bin/gff/main.py
	python $^ -l 250 -ml 15 -f 4 -mf 4 -loc "data/families/WGD2/upstream/" --head

retrieve_CDS: bin/gff/ntseq.py
	python $^ -f 4 --header 
