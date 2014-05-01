# Make test
# .PHONY: all

all: WGD2ANC00182.newick


WGD2ANC%.newick: WGD2ANC%CDS.phyl_phyml_tree.txt
	../../bin/scripts/editnewick.py $< $@

WGD2ANC%CDS.phyl_phyml_tree.txt: WGD2ANC%CDS.phyl
	phyml -i $<

WGD2ANC%CDS.phyl: WGD2ANC%CDS.e.fasta
	perl ../../bin/scripts/ConvertFastatoPhylip.pl WGD2ANC00208CDS.e.fasta WGD2ANC00208CDS.phyl


