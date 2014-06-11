#!/usr/bin/env Rscript
# highnonmatch.R
# Find non matching names between highly expressed and ribosomal genes

### Libraries ###
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(argparse))


### Functions ###
high.non.match <- function(expression_file, ribo.seqs, quantile) {
# Extract names of ribosomal genes
expression_file = read.table(expression_file, h=TRUE)
ribo.seqs = read.fasta(ribo.seqs)
ribo.names = attr(ribo.seqs, "name")

# Compute quantile 
print(quantile(expression_file$xp, probs=seq(0.9,1,0.01)))
quant = quantile(expression_file$xp, probs=quantile)

threshold = as.numeric(quant)  # Take 98% quantile

high_exp = subset(expression_file, xp >= threshold)

# Return list of non ribosomal highly expressed genes
matching = high_exp[!high_exp$gene_id %in% as.factor(ribo.names), 1]
}


### Command-line parser ###
parser <- ArgumentParser(description="Script to extract the names highly expressed non ribosomal genes. Extract top expressed genes according to threshold, from expression file not matching sequences file.") 
parser$add_argument("expfile", help="RNAseq expression tabfile, with a column \'xp\' and a column \'gene_id\' with gene name")
parser$add_argument("sequences", help="Fasta file with the name of each gene as header")
parser$add_argument("-t", "--thre", type="double",  default=0.98, help="Consider only genes over this percentile. Example: if t = 0.98, consider only top 2 percent highly expressed genes. (default %(default)s)")
parser$add_argument("-o", "--output", default="highnonmatch.txt", help="name of the output file (default %(default)s)")
args <- parser$parse_args()


### Data loading ###
exp = read.table(args$expfile, h=TRUE)  # Expression data
ribo.seqs = read.fasta(args$sequences)  # Ribosomal upstream sequences


### Computing ###
matching = high.non.match(exp, ribo.seqs, args$thre)
print(paste("Writing file ", args$output))
write.table(matching, args$output, quote=FALSE, row.names=FALSE, col.names=FALSE)
print("Done.")