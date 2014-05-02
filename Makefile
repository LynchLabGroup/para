# Pipeline make file
# with all the steps

# percent_subdirs := 1 # Allow % to match multiple directories
SUBDIRS = $(shell find results/ -type d -name 'WGD2ANC0000*')
FILES = $(patsubst results/WGD2ANC%, results/WGD2ANC%/WGD2ANC\%, $(SUBDIRS))
FAM = $(subst results/, ,$(SUBDIRS))
.PHONY: all


all: 
	@echo $(SUBDIRS)
	@echo $(FILES)
	@echo $(FAM)


# WGD2ANC%.newick: WGD2ANC%CDS.phyl_phyml_tree.txt
# 	../../bin/scripts/editnewick.py $< $@

# WGD2ANC%CDS.phyl_phyml_tree.txt: WGD2ANC%CDS.phyl
# 	phyml -i $<

# WGD2ANC%CDS.phyl: WGD2ANC%CDS.e.fasta
# 	perl ../../bin/scripts/ConvertFastatoPhylip.pl WGD2ANC00208CDS.e.fasta WGD2ANC00208CDS.phyl


