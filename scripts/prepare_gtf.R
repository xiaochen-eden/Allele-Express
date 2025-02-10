library(tximport)
library(readr)
library(GenomicFeatures)
library(AnnotationDbi)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("请提供3个参数：GTF文件路径，转录本长度文件路径，输出文件路径")
}

gtf_file <- args[1]
transcript_length_file <- args[2]
output_file <- args[3]

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, columns = c("GENEID"), keytype = "TXNAME")
transcript_lengths <- read.table(transcript_length_file, header = FALSE, col.names = c("TXNAME", "Length"))
tx_info <- merge(tx2gene, transcript_lengths, by = "TXNAME")
write.table(tx_info, output_file, row.names = FALSE, sep = "\t", quote = FALSE)

