#!/usr/bin/env Rscript

args=(commandArgs(TRUE))
regenie_setinclusionfile=as.character(args[1])
gene_chunk_num=as.numeric(args[2])
gene_chunk_total=as.numeric(args[3])
tissues_to_keep=as.character(args[4])
regenie_setinclusionfile_out=as.character(args[5])

library(data.table)
.libPaths(c("rpackages4_1_3",.libPaths()))

group <- fread(regenie_setinclusionfile, stringsAsFactors = F, data.table=F, header=F, sep='\t')
genes <- unique(gsub("__.*", "", group$V1))

## Find chunk genes to run
chunk_size <- ceiling(length(genes)/gene_chunk_total)
splitz <- base::split(c(1:length(genes)), ceiling(seq_along(c(1:length(genes)))/chunk_size))
genes_chunk <- genes[splitz[[gene_chunk_num]]]

## Filter to chunk genes
group <- as.data.frame(group[which(gsub("__.*", "", group$V1) %in% genes_chunk), ], stringsAsFactors=FALSE)
colnames(group) <- "V1"

## Filter to needed tissues
tissues_to_keep <- strsplit(tissues_to_keep, split=",")
tissues_to_keep_vec <- NULL
for(i in c(1:length(tissues_to_keep[[1]]))){
  tissues_to_keep_vec <- c(tissues_to_keep_vec, tissues_to_keep[[1]][i])
}
group <- as.data.frame(group[which(gsub(".*__", "", group$V1) %in% tissues_to_keep_vec), ], stringsAsFactors=FALSE)
colnames(group) <- "V1"

## Order by gene
group <- group[order(group$V1), ]

## Save result
write.table(group, file=regenie_setinclusionfile_out, col.names=F, row.names=F, quote=F, sep='\t')
