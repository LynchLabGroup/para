# highnonmatch.R
# Find non matching names between highly expressed and ribosomal genes

library(seqinr)

# Data loading
exp = read.table("data/tetraurelia_RNAseq_xp.tab", h=T) # Expression data
ribo.seqs = read.fasta("../ribo.upstream.fa") # Ribosomal upstream sequences

# Extract names of ribosomal genes
ribo.names = attr(ribo.seqs, "name")

# Compute quantile 
quant = quantile(exp$xp, probs=seq(0.9,1,0.01))

threshold = quant[9] # Take 98% quantile

high_exp = subset(exp, xp >= threshold)

# Return list of non ribosomal highly expressed genes
matching = high_exp[!high_exp %in% as.factor(ribo.names)]

write.table(matching, "results/highnonmatch.txt")