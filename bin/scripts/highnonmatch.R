# highnonmatch.R
# Find non matching names between highly expressed and ribosomal genes

library(seqinr)

# Data loading
exp = read.table("data/tetraurelia_RNAseq_xp.tab", h=TRUE) # Expression data
ribo.seqs = read.fasta("../ribo.upstream.fa") # Ribosomal upstream sequences

high.non.match <- function(expression_file, ribo.seqs) {
# Extract names of ribosomal genes
expression_file = read.table(expression_file, h=TRUE)
ribo.seqs = read.fasta(ribo.seqs)
ribo.names = attr(ribo.seqs, "name")

# Compute quantile 
print(quantile(expression_file$xp, probs=seq(0.9,1,0.01)))
quant = quantile(expression_file$xp, probs=seq(0.9,1,0.01))

threshold = as.numeric(quant[9]) # Take 98% quantile

high_exp = subset(expression_file, xp >= threshold)

# Return list of non ribosomal highly expressed genes
matching = high_exp[!high_exp[,1] %in% as.factor(ribo.names), 1]

write.table(matching, "results/highnonmatch.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
}