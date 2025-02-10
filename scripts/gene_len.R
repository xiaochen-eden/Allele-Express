#!/usr/bin/env Rscript

# ---------------------------------------------------------
# Script Name: gene_len.R
# Author: chenxiao
# Date: 2024/07/16
# Description: This script performs norm analysis on MAS sequencing data.
# ---------------------------------------------------------

library(parallel)
cl <- makeCluster(0.75*detectCores())
library(GenomicFeatures)

args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[1]
gene_quant_file <- args[2]
output_file <- args[3]

gene_ids <- read.table(gene_quant_file, header = FALSE, sep = "\t")[, 1]

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
exons_gene <- exonsBy(txdb, by = "gene")

exons_gene_lens <- lapply(exons_gene, function(x) { sum(width(reduce(x))) })
geneid_efflen <- data.frame(geneid = names(exons_gene_lens), efflen = as.numeric(exons_gene_lens))

# Filter to ensure only desired gene IDs are present
geneid_efflen <- geneid_efflen[geneid_efflen$geneid %in% gene_ids, ]

write.table(geneid_efflen, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

