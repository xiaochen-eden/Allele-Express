#!/usr/bin/env Rscript

# ---------------------------------------------------------
# Script Name: norml.R
# Author: chenxiao
# Date: 2024/07/16
# Description: This script performs norm analysis on MAS sequencing data.
# ---------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

a1 <- fread(input_file, header = TRUE, data.table = FALSE)
a1 <- a1[!grepl("__ambiguous|__no_feature|__not_aligned|geneid", a1$Geneid), ]
print(a1)

counts <- a1[, 7:ncol(a1)]
counts_df <- data.frame(counts = counts)
rownames(counts_df) <- a1$Geneid

geneid_efflen <- subset(a1, select = c("Geneid", "Length"))
colnames(geneid_efflen) <- c("geneid", "efflen")
efflen <- geneid_efflen[match(rownames(counts_df), geneid_efflen$geneid), "efflen"]
print(efflen)

efflen <- efflen[efflen != "" & efflen != "efflen"] 
efflen_numeric <- as.numeric(efflen)  

### 计算 TPM
#TPM (Transcripts Per Kilobase Million)  每千个碱基的转录每百万映射读取的Transcripts
counts2TPM <- function(count, efflength){
  RPK <- count/(efflength/1000)       #每千碱基reads (“per million” scaling factor) 长度标准化
  PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
  RPK/PMSC_rpk                    
}  
tpm <- as.data.frame(apply(counts_df, 2, counts2TPM,efflen_numeric))
write.table(tpm, output_file, sep = "\t",col.names = FALSE, quote = FALSE)

