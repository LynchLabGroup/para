# Make test
.PHONY: all clean

all : WGD2ANC00182.newick WGD2ANC00208.newick

WGD2ANC%.newick: WGD2ANC%CDS.phyl_phyml_tree.txt
	python ../../bin/scripts/editnewick.py $< $@

WGD2ANC%CDS.phyl_phyml_tree.txt: WGD2ANC%CDS.phyl
	phyml -i $<

WGD2ANC%CDS.phyl: WGD2ANC%CDS.e.fasta
	perl ../../bin/scripts/ConvertFastatoPhylip.pl $< $@

clean:
	rm *.txt
	rm *.phyl
	rm *.newick
