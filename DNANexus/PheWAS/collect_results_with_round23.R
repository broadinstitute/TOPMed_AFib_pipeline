#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
keyfile=as.character(args[1])
overvfile=as.character(args[2])
chunk_size=as.numeric(args[3])
chunk_num=as.numeric(args[4])

.libPaths(c("rpackages4_1_3",.libPaths()))

#git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#git pull --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#source("/medpop/afib/sjurgens/Rscripts/association_source_v2.R")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("UKBB_200KWES_CVD/Cauchy_test.R")
install.packages("R.utils")
library(data.table)
library(dplyr)
library(tidyr)

#############################
# Key file and overview files
#############################

## Read in key file to find the phenotype
key <- fread(keyfile, stringsAsFactors=F, data.table=F)
key <- key[key$Ancestry=="ALL", ]
key <- key[key$N_cases >= 50 & key$N_controls >=50, ]
dim(key)

overv <- fread(overvfile, stringsAsFactors = F, data.table=F)
overv <- overv[overv$included_in_mgb_phewas=="yes", ]

## Find chunk to run
#splitz <- split(c(1:nrow(key)), ceiling(seq_along(c(1:nrow(key)))/chunk_size))
#key <- key[splitz[[chunk_num]], ]
#key <- key[chunk_num, ]
i <- chunk_num

##### Check outfiles exist #####
#foreach(i=c(1:nrow(key)), .inorder=FALSE) %dopar% {
#for(i in c(1:nrow(key))){ 
inter <- inter2 <- NULL
num <- key[i, 'Phecode']
n.cases <- key[i, 'N_cases']
n.controls <- key[i, 'N_controls']
phenoname <- key[i, 'Name']
category <- overv[overv$meaning==phenoname, 'category']
cat('\n\n\nBusy with phenotype', num, 'which is', phenoname, 'and task', i, 'out of 535 tasks...\n\n')

######################
# Result files
######################

cat('\tchecking if results files are present...\n\n')
#MAF<0.1% files
#for(chr in c(1:22)){
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_lowmem.RData")))
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_highmem.RData")))
#    if(chr %in% c(1:20, 22)){
#        try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_highhighmem.RData")))
#    }
#    if(chr %in% c(1:8, 10:13, 15:17, 19:20)){
#        try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_veryhighmem.RData")))
#    }
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/", num, "_results_chr", chr, "_maf0.001.RData")))
#}
files <- paste0("dx download")
for(chr in c(1:22)){
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_lowmem.RData"))
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_highmem.RData"))
    if(chr %in% c(1:20, 22)){
        files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_highhighmem.RData"))
    }
    if(chr %in% c(1:8, 10:13, 15:17, 19:20)){
        files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.001_round2_veryhighmem.RData"))
    }
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/", num, "_results_chr", chr, "_maf0.001.RData"))
}
#files1_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/lowmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.001_round2_lowmem.RData")
#files1_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/highmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.001_round2_highmem.RData")
#files1_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/highhighmem/chr", c(1:20, 22), '/', num, "_results_chr", c(1:20, 22), "_maf0.001_round2_highhighmem.RData")
#files1_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/veryhighmem/chr", c(1:8, 10:13, 15:17, 19:20), '/', num, "_results_chr", c(1:8, 10:13, 15:17, 19:20), "_maf0.001_round2_veryhighmem.RData")
#files1_5 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/", num, "_results_chr", c(1:22), "_maf0.001.RData")
files1_1 <- paste0(num, "_results_chr", c(1:22), "_maf0.001_round2_lowmem.RData")
files1_2 <- paste0(num, "_results_chr", c(1:22), "_maf0.001_round2_highmem.RData")
files1_3 <- paste0(num, "_results_chr", c(1:20, 22), "_maf0.001_round2_highhighmem.RData")
files1_4 <- paste0(num, "_results_chr", c(1:8, 10:13, 15:17, 19:20), "_maf0.001_round2_veryhighmem.RData")
files1_5 <- paste0(num, "_results_chr", c(1:22), "_maf0.001.RData")
maf0.001_files <- c(files1_1, files1_2, files1_3, files1_4, files1_5)
maf0.001_nfilesets <- 5

#MAF<0.001% files
#for(chr in c(1:22)){
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_lowmem.RData")))
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_highmem.RData")))
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_highhighmem.RData")))
#    if(chr %in% c(1:17, 19:20)){
#        try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round2/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_veryhighmem.RData")))
#    }
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_lowmem.RData")))
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_highmem.RData")))
#    if(chr %in% c(1:20, 22)){
#        try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_highhighmem.RData")))
#    }
#    if(chr %in% c(1:6, 8, 12:13, 15:16, 19)){
#        try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_veryhighmem.RData")))
#    }
#    # Old MAF<0.001% run with some errors
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/", num, "_results_chr", chr, "_maf0.00001.RData")))
#}
for(chr in c(1:22)){
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_lowmem.RData"))
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_highmem.RData"))
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_highhighmem.RData"))
    if(chr %in% c(1:17, 19:20)){
        files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round2/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round2_veryhighmem.RData"))
    }
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_lowmem.RData"))
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_highmem.RData"))
    if(chr %in% c(1:20, 22)){
        files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_highhighmem.RData"))
    }
    if(chr %in% c(1:6, 8, 12:13, 15:16, 19)){
        files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.00001_round3_veryhighmem.RData"))
    }
    # Old MAF<0.001% run with some errors
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/", num, "_results_chr", chr, "_maf0.00001.RData"))
}
#files2_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/lowmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.00001_round2_lowmem.RData")
#files2_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/highmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.00001_round2_highmem.RData")
#files2_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/highhighmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.00001_round2_highhighmem.RData")
#files2_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/veryhighmem/chr", c(1:17, 19:20), '/', num, "_results_chr", c(1:17, 19:20), "_maf0.00001_round2_veryhighmem.RData")
#files2_5 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/lowmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.00001_round3_lowmem.RData")
#files2_6 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/highmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.00001_round3_highmem.RData")
#files2_7 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/highhighmem/chr", c(1:20, 22), '/', num, "_results_chr", c(1:20, 22), "_maf0.00001_round3_highhighmem.RData")
#files2_8 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/veryhighmem/chr", c(1:6, 8, 12:13, 15:16, 19), '/', num, "_results_chr", c(1:6, 8, 12:13, 15:16, 19), "_maf0.00001_round3_veryhighmem.RData")
files2_1 <- paste0(num, "_results_chr", c(1:22), "_maf0.00001_round2_lowmem.RData")
files2_2 <- paste0(num, "_results_chr", c(1:22), "_maf0.00001_round2_highmem.RData")
files2_3 <- paste0(num, "_results_chr", c(1:22), "_maf0.00001_round2_highhighmem.RData")
files2_4 <- paste0(num, "_results_chr", c(1:17, 19:20), "_maf0.00001_round2_veryhighmem.RData")
files2_5 <- paste0(num, "_results_chr", c(1:22), "_maf0.00001_round3_lowmem.RData")
files2_6 <- paste0(num, "_results_chr", c(1:22), "_maf0.00001_round3_highmem.RData")
files2_7 <- paste0(num, "_results_chr", c(1:20, 22), "_maf0.00001_round3_highhighmem.RData")
files2_8 <- paste0(num, "_results_chr", c(1:6, 8, 12:13, 15:16, 19), "_maf0.00001_round3_veryhighmem.RData")
maf0.00001_files <- c(files2_1, files2_2, files2_3, files2_4, files2_5, files2_6, files2_7, files2_8)
maf0.00001_nfilesets <- 8
# Old MAF<0.001% run with some errors
#files2_9 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/", num, "_results_chr", c(1:22), "_maf0.00001.RData")
files2_9 <- paste0(num, "_results_chr", c(1:22), "_maf0.00001.RData")
maf0.00001_old_files <- c(files2_9)
maf0.00001_old_nfilesets <- 1

#MAF<1% files
#for(chr in c(1:22)){
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_lowmem.RData")))
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_highmem.RData")))
#    try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_highhighmem.RData")))
#    if(chr %in% c(1:20, 22)){
#        try(system(paste0("dx download exome-seq:/sjj/projects/phewas/v1/results/association/round3/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_veryhighmem.RData")))
#    }
#}
for(chr in c(1:22)){
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/lowmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_lowmem.RData"))
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/highmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_highmem.RData"))
    files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/highhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_highhighmem.RData"))
    if(chr %in% c(1:20, 22)){
        files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/round3/veryhighmem/chr", chr, '/', num, "_results_chr", chr, "_maf0.01_round3_veryhighmem.RData"))
    }
}
#files3_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/lowmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.01_round3_lowmem.RData")
#files3_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/highmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.01_round3_highmem.RData")
#files3_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/highhighmem/chr", c(1:22), '/', num, "_results_chr", c(1:22), "_maf0.01_round3_highhighmem.RData")
#files3_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round3/veryhighmem/chr", c(1:20, 22), '/', num, "_results_chr", c(1:20, 22), "_maf0.01_round3_veryhighmem.RData")
files3_1 <- paste0(num, "_results_chr", c(1:22), "_maf0.01_round3_lowmem.RData")
files3_2 <- paste0(num, "_results_chr", c(1:22), "_maf0.01_round3_highmem.RData")
files3_3 <- paste0(num, "_results_chr", c(1:22), "_maf0.01_round3_highhighmem.RData")
files3_4 <- paste0(num, "_results_chr", c(1:20, 22), "_maf0.01_round3_veryhighmem.RData")
maf0.01_files <- c(files3_1, files3_2, files3_3, files3_4)
maf0.01_nfilesets <- 4

# Dx downloading the list of files...
try(system(paste0(files)))

if(!(all(file.exists(maf0.001_files)) & all(file.exists(maf0.00001_files)) & all(file.exists(maf0.01_files)))){
    cat(paste0(which(!file.exists(c(maf0.001_files, maf0.00001_files, maf0.01_files))), "\n"))
    cat("\n\n\n\n\nWARNING: not all required files found!!! Stopping.\n\n\n\n")
}else{
    cat('\tall results files found. Reading in and merging...\n\n')
    cat('\t\tMAF<0.1%...\n\n')
    
    #MAF<0.1%
    inter <- summarydata(files=maf0.001_files, chrs=c(c(1:22), c(1:22), c(1:20, 22), c(1:8, 10:13, 15:17, 19:20), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- rownames(inter)
    inter$gene <- sub("_", "__", inter$gene)
    inter$gene <- sub("_gnomAD_POPMAX0.001", "_POPMAX0.001", inter$gene)
    
    cat('\t\tMAF<0.001%...\n\n')
    #MAF<0.001% 
    inter2 <- summarydata(files=maf0.00001_files, chrs=c(c(1:22), c(1:22), c(1:22), c(1:17, 19:20), c(1:22), c(1:22), c(1:20, 22), c(1:6, 8, 12:13, 15:16, 19)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter2 <- inter2[inter2$n.sample.alt>=10, ]
    inter2$category <- category
    inter2$gene <- rownames(inter2)
    inter2$gene <- sub("_", "__", inter2$gene)
    inter2$gene <- sub("_gnomAD_POPMAX0.00001", "_POPMAX0.00001", inter2$gene)
    #Old (with some errors)
    if(any(file.exists(maf0.00001_old_files))){
      inter2_2 <- summarydata(files=maf0.00001_files, chrs=rep(c(1:22),maf0.00001_nfilesets), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
      inter2_2 <- inter2_2[inter2_2$n.sample.alt>=10, ]
      inter2_2$category <- category
      inter2_2$gene <- rownames(inter2_2)
      inter2_2$gene <- sub("_", "__", inter2_2$gene)
      inter2_2$gene <- sub("_gnomAD_POPMAX0.00001", "_POPMAX0.00001", inter2_2$gene)
      inter2 <- rbind(inter2, inter2_2)
      rm <- which(duplicated(inter2$gene))
      if(length(rm)>0){inter2 <- inter2[-rm, ]}
    }
    
    cat('\t\tMAF<1%...\n\n')
    #MAF<1%
    inter3 <- summarydata(files=maf0.01_files, chrs=c(c(1:22), c(1:22), c(1:22), c(1:20, 22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter3 <- inter3[inter3$n.sample.alt>=10, ]
    inter3$category <- category
    inter3$gene <- rownames(inter3)
    inter3$gene <- sub("_", "__", inter3$gene)
    inter3$gene <- sub("_gnomAD_POPMAX0.01", "_POPMAX0.01", inter3$gene)

    #Combine
    inter <- rbind(inter, inter2, inter3)
    inter$mask <- gsub(".*__", "", inter$gene)
    inter$gene <- gsub("__.*", "", inter$gene)
    inter$cases <- n.cases
    inter$controls <- n.controls
                           
    rawassoc_res <- inter
    write.table(rawassoc_res, file=paste0('../summary_results_phewas_all_tests_phecode', num, '.tsv'), col.names=F, row.names=F, quote=F, sep='\t', append=T)

    # Create Cauchy file; run first only for truly rare variant masks
    cat('\trunning Cauchy combinations for MAF<0.1% and 0.001% masks...\n\n')
    inter <- inter[inter$n.sample.alt>=20, ]
    inter1 <- inter[inter$mask=="hclof_noflag_POPMAX0.001", c("phenotype", "cases", "controls", "gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter1)[c(5:9)] <- paste0(colnames(inter1)[c(5:9)], "__hclof_noflag_POPMAX0.001")
    inter2 <- inter[inter$mask=="hclof_noflag_missense0.8_POPMAX0.001", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter2)[c(2:6)] <- paste0(colnames(inter2)[c(2:6)], "__hclof_noflag_missense0.8_POPMAX0.001")
    inter3 <- inter[inter$mask=="hclof_noflag_missense0.5_POPMAX0.001", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter3)[c(2:6)] <- paste0(colnames(inter3)[c(2:6)], "__hclof_noflag_missense0.5_POPMAX0.001")
    inter4 <- inter[inter$mask=="hclof_noflag_missense0.5_POPMAX0.00001", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter4)[c(2:6)] <- paste0(colnames(inter4)[c(2:6)], "__hclof_noflag_missense0.5_POPMAX0.00001")
    inter5 <- inter[inter$mask=="missense0.5_POPMAX0.00001", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter5)[c(2:6)] <- paste0(colnames(inter5)[c(2:6)], "__missense0.5_POPMAX0.00001")
    inter6 <- inter[inter$mask=="missense0.2_POPMAX0.00001", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter6)[c(2:6)] <- paste0(colnames(inter6)[c(2:6)], "__missense0.2_POPMAX0.00001")

    inter1 <- merge(inter1, inter2, by="gene", all=T)
    inter1 <- merge(inter1, inter3, by="gene", all=T)
    inter1 <- merge(inter1, inter4, by="gene", all=T)
    inter1 <- merge(inter1, inter5, by="gene", all=T)
    inter1 <- merge(inter1, inter6, by="gene", all=T)
    colz <- which(grepl("SPA.pval", colnames(inter1)))
    inter1$P_cauchy  <- apply(inter1[,colz], 1, CCT)

    # Including MAF<1%
    cat('\trunning Cauchy combinations including also the MAF<1% masks...\n\n')
    inter2 <- inter[inter$mask=="hclof_noflag_POPMAX0.01", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter2)[c(2:6)] <- paste0(colnames(inter2)[c(2:6)], "__hclof_noflag_POPMAX0.01")
    inter3 <- inter[inter$mask=="hclof_noflag_missense0.8_POPMAX0.01", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter3)[c(2:6)] <- paste0(colnames(inter3)[c(2:6)], "__hclof_noflag_missense0.8_POPMAX0.01")
    inter4 <- inter[inter$mask=="hclof_noflag_missense0.5_POPMAX0.01", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter4)[c(2:6)] <- paste0(colnames(inter4)[c(2:6)], "__hclof_noflag_missense0.5_POPMAX0.01")

    inter1 <- merge(inter1, inter2, by="gene", all=T)
    inter1 <- merge(inter1, inter3, by="gene", all=T)
    inter1 <- merge(inter1, inter4, by="gene", all=T)
    colz <- which(grepl("SPA.pval", colnames(inter1)))
    inter1$P_cauchy_v2  <- apply(inter1[,colz], 1, CCT)

    inter1$phenotype <- phenoname
    inter1$cases <- n.cases
    inter1$controls <- n.controls
    write.table(inter1, file=paste0('../summary_results_phewas_cauchy_phecode', num, '.tsv'), col.names=F, row.names=F, quote=F, sep='\t', append=T)

  #############################
  # Build phenotype files 
  #############################

  cat('\tcreating phenotype and helper files for analysis...\n\n')
  # Example phenotype file from the same phenotype freeze
  exdat0 <- fread("/mnt/project/sjj/projects/phewas/v1/data/pheno/Hypertrophic_cardiomyopathy.tab.tsv.gz",header=T,data.table=F)
  linker <- fread("/mnt/project/sjj/projects/phewas/v1/data/pheno/ukb_app17488_app7089_link.csv",header=T,data.table=F,sep=",")
  exdat1 <- merge(exdat0,linker,by.x="sample_id",by.y="app7089")
  exdat1 <- exdat1[,c('sample_id', 'enroll_age', 'has_died', 'app17488')]
  exdat1$enroll_age_2 <- (exdat1$enroll_age)^2
  head(exdat1)
  dim(exdat1)

  # Plink sample files
  ##### 450k is the sampleQCd file
  u450k <- fread('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c1_genotype_variant_sample_QCed.psam',
              stringsAsFactors=F, data.table=F, header=F)
  u450k$u450k <- 1
  u450k <- u450k[,c("V1", "u450k")]
  length(which(u450k$V1 %in% exdat1$sample_id))
  length(which(u450k$V1 %in% exdat1$app17488))
  ### app17488 is genetic data, as it should

  ##### 200k and 50k files are raw fam files from previous freezes
  u200k <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_200k_tranche.fam',
              stringsAsFactors=F, data.table=F, header=F)
  u200k$u200k <- 1
  u200k <- u200k[,c("V1", "u200k")]
  length(which(u200k$V1 %in% exdat1$sample_id))
  length(which(u200k$V1 %in% exdat1$app17488))
  ### 200k is on app17488 as well

  u50k <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_50k_tranche.fam',
              stringsAsFactors=F, data.table=F, header=F)
  u50k$u50k <- 1
  u50k <- u50k[,c("V1", "u50k")]
  length(which(u50k$V1 %in% exdat1$sample_id))
  length(which(u50k$V1 %in% exdat1$app17488))
  ### 50k is on other app

  phen0 <- merge(exdat1, u450k, by.x='app17488', by.y='V1', all=F)
  phen0 <- merge(phen0, u200k, by.x='app17488', by.y='V1', all.x=T, all.y=F)
  phen0[is.na(phen0$u200k), 'u200k'] <- 0
  phen0 <- merge(phen0, u50k, by.x='sample_id', by.y='V1', all.x=T, all.y=F)
  phen0[is.na(phen0$u50k), 'u50k'] <- 0
  head(phen0)
  table(phen0$u450k)
  table(phen0$u200k)
  table(phen0$u50k)

  # Batch variable
  phen0$batch <- "u450k"
  phen0[phen0$u200k==1, 'batch'] <- "u200k"
  phen0[phen0$u50k==1, 'batch'] <- "u50k"
  table(phen0$batch)

  # Admixture, PCA and sex check results
  adm <- fread('/mnt/project/exome_450k_plink/ADMIXTURE/ADMIXTURE_results_UKBB450k.txt',
               stringsAsFactors=F, data.table=F)
  adm <- adm[,c("ID", "POP_tight")]
  adm$ID <- gsub("_.*", "", adm$ID)
  pca <- fread('/mnt/project/exome_450k_plink/PCA/PCs_round1.txt',
            stringsAsFactors=F, data.table=F)
  pca$sample_id <- gsub("_.*", "", pca$sample_id)
  sexc <- fread('/mnt/project/exome_450k_plink/sampleQC/sex_check/sex_mismatch_report_complete.tsv',
             stringsAsFactors=F, data.table=F)
  sexc <- sexc[,c("IID", "reported_gender")]
  colnames(sexc)[2] <- "Gender"

  phen0 <- merge(phen0, sexc, by.x='app17488', by.y='IID')
  phen0 <- merge(phen0, adm, by.x="app17488", "ID")
  phen0 <- merge(phen0, pca, by.x="app17488", by.y="sample_id")
  dim(phen0)
  head(phen0)

  # Raw phenotype case file; James' file follows other ID
  raw <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_phecode_affectedsamples_raw.csv.gz', 
             stringsAsFactors=F, data.table=F)
  class(raw$jd_code) <- "numeric"
  #length(which(exdat1$app17488 %in% raw$sample_id))
  #length(which(exdat1$sample_id %in% raw$sample_id))

  # Some additional phenotypes were manually curated by me; will add if necessary
  if(chunk_num %in% c((550-14):(560-14))){
    extra_phens <- fread('/mnt/project/sjj/projects/phewas/v1/data/pheno/ukbb_missing_phecodes_additionally_curated.tsv', stringsAsFactors=F, data.table=F)  
    phen0 <- merge(phen0, extra_phens, by.x="app17488", by.y="ID", all.x=T)
  }
  
  #system('dx download exome-seq:/exome_450k_plink/PCA/ukb23156_KING_GRM_sparseto3rddegree_scaledby2.RData')
  #mat <- get(load('ukb23156_KING_GRM_sparseto3rddegree_scaledby2.RData'))
  #colnames(mat) <- gsub("_.*", "", colnames(mat))
  #rownames(mat) <- gsub("_.*", "", rownames(mat))
  #mat[1:10,1:10]
  #save(mat, file='ukb23156_KING_GRM_sparseto3rddegree_scaledby2.RData')

  jd_code <- overv[overv$meaning==phenoname, 'jd_ukbb_code']
  raw <- raw[raw$jd_code ==jd_code, ]

  # Define the disease variable based on the disease files
  phen0$disease <- 0
  if(chunk_num %in% c((550-14):(560-14))){
    phen0$disease <- phen0[,paste0('phecode_', jd_code)]
  }else{
    # James' IDs are coded for the other application
    phen0[phen0$sample_id %in% raw$sample_id, 'disease'] <- 1
  }
  gender <- overv[overv$meaning==phenoname, 'Gender']
  rm <- which(is.na(phen0$app17488))
  if(length(rm)>0){phen0 <- phen0[ -rm,]}

  #impute missing to 0
  phen0[is.na(phen0[,'disease']), 2] <- 0
        
  if(gender == "male"){
        phen0 <- phen0[phen0$Gender=="Male" & !is.na(phen0$Gender), ]
  }else if(gender == "female"){
        phen0 <- phen0[phen0$Gender=="Female" & !is.na(phen0$Gender), ]
  }
  
  #nullmodel file to extract covariates; categorial variables need to be converted 
  nullmod <- get(load(paste0('/mnt/project/sjj/projects/phewas/v1/nullmodels/', jd_code, '_nullmodel.RData')))
  fixef <- colnames(nullmod$betaCov)[c(2:ncol(nullmod$betaCov))]
  catCovarList <- NULL
  if(any(grepl("Gender", fixef))){
      fixef[which(grepl("Gender", fixef))] <- "Gender"
      catCovarList <- c(catCovarList, "Gender")
  }
  if(any(grepl("batch", fixef))){
      fixef[which(grepl("batch", fixef))] <- "batch"
      catCovarList <- c(catCovarList, "batch")
  }
  fixef <- unique(fixef)
  cat("\t\tcovariates of interest:\n")
  print(fixef)
  cat("\n\n")
    
  # filter to samples used in GENESIS nullmod; this is important for variant filtering in PLINK
  nullmod_samples <- nullmod$sample.id
  phen0 <- phen0[phen0$app17488 %in% nullmod_samples, ]
  write.table(cbind(phen0$app17488,phen0$app17488), file=paste0(num, '__sampleIDs.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  phen0 <- phen0[ ,c(1, 1, c(2:(ncol(phen0))))]
  colnames(phen0)[c(1,2)] <- c("FID", "IID")
  write.table(phen0, file=paste0(num, '__regenie_phenofile.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
  write.table(phen0[phen0$POP_tight=="EUR", c("FID", "IID")], file=paste0(num, '__sampleIDs_EUR.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  cat(paste0("\t\ttotal sample size is: ", nrow(phen0), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1, ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0, ]), "\n\n\n"))

  # For running Firth's we need to also define unrelated samples for running in REGENIE...
  system('dx download exome-seq:/exome_450k_plink/PCA/ukbb_450k_unrelatedsamples.tsv')
  unrel <- fread('ukbb_450k_unrelatedsamples.tsv', stringsAsFactors = F, data.table=F)
  colnames(unrel) <- "sample.id"
  phen_unrel <- phen0[phen0$IID %in% unrel$sample.id, ]
  unrel_samples <- cbind(phen_unrel$FID,phen_unrel$IID)
  colnames(unrel_samples) <- c("FID", "IID")
  write.table(unrel_samples, file=paste0(num, '__sampleIDs_unrel.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
  firth.n.cases <- nrow(phen_unrel[phen_unrel$disease==1, ])
  firth.n.controls <- nrow(phen_unrel[phen_unrel$disease==0, ])
  cat(paste0("\t\tunrelated sample size is: ", nrow(phen_unrel), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls, "\n\n\n"))

  ##########################################
  ## Run PLINK -> REGENIE Firth pipeline
  ##########################################

  ###### Identify tests to rerun
  cat('\tRunning PLINK extracting and REGENIE effect size estimation...\n\n')
  ### identify based on whetehr SPA was applied (if SPA was applied, then the nominal P<0.05); also rerun P<0.1 with beta>0
  assocs <- rawassoc_res[(rawassoc_res$SPA.converged & !is.na(rawassoc_res$SPA.converged)) | (rawassoc_res$SPA.pval<0.1 & rawassoc_res$Est>0), ]
  gene_masks <- paste0(assocs$gene, "__", assocs$mask)
  
  regenie_res_tot <- NULL
  regenie_res_tot_nonEUR <- NULL

  ######### MAF<0.1% masks
  #Create set files for regenie; filter PLINK files to needed variants only!
  regenie_setfile <- NULL
  regenie_annotationfile <- NULL
  cat("\t\textracting variant data from PLINK files for MAF<0.1% threshold ...\n")

  for(chr in c(1:22)){
        cat("\t\t\tbusy with chromosome ", chr, "...\n")
        plink_path <- './plink2'
        regenie_path <- './regenie_v3.2.2.gz_x86_64_Linux_mkl'
        plinkfile <- paste0('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c', chr, '_genotype_variant_sample_QCed')
        plinkfile_type <- "pfile"
        max_maf=0.001
        max_mac='100000000'
        #carz <- NULL
        #phen1 <- phen0

        ##Use the grouping files used for the analyses
        files1_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_lowmem.RData")
        files1_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highmem.RData")
        files1_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highhighmem.RData")
        files1_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_veryhighmem.RData")
        files1_5 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_popmax0.001.RData")

        f1 <- get(load(files1_1))
        f2 <- get(load(files1_2)) 
        f3 <- get(load(files1_3)) 
        f4 <- get(load(files1_4)) 
        f5 <- get(load(files1_5))
        class(f1$pos) <- class(f2$pos) <- class(f3$pos) <- class(f4$pos) <- class(f5$pos) <- "numeric"
        class(f1$chr) <- class(f2$chr) <- class(f3$chr) <- class(f4$chr) <- class(f5$chr) <- "numeric"
        group <- bind_rows(f1, f2, f3, f4, f5)

        group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
        group$group_id <- sub("_", "__", group$group_id)
        group$group_id <- sub("_gnomAD_POPMAX0.001", "_POPMAX0.001", group$group_id)

        ## Filter to the masks we want to test
        group <- group[group$group_id %in% gene_masks, ]
        #group <- group[,c("varid", "alt", "group_id")]
        group$varid <- paste0("chr", group$varid)

        ## Make part of REGENIE grouping files; collect groupings; will use later on
        group$pseudo_annot <- "REGENIE"
        regenie_annotationfile <- rbind(regenie_annotationfile, group[,c("varid", "group_id", "pseudo_annot")])
        set_inter <- NULL
        for(gr in unique(group$group_id)){
                gr_inter <- group[group$group_id==gr, ]
                chromosome <- gr_inter[1,'chr']
                position <- gr_inter[1,'pos']
                collapse <- paste(gr_inter$varid, collapse=",")
                set_inter <- rbind(set_inter, c(gr, chromosome, position, collapse))
        }
        regenie_setfile <- rbind(regenie_setfile, set_inter)
        
        ## Use PLINK2 to filter to a PLINK file that has the needed variants and samples used in the discovery analysis
        write.table(group$varid, file=paste0(num, '__varz_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
        try(system(paste0(plink_path, ' ',
                '--', plinkfile_type, '  ', plinkfile, '  ',
                '--max-maf ', max_maf, '  --max-mac ', max_mac, '  ',
                '--keep  ', num, '__sampleIDs.tsv  ',
                '--extract  ', num, '__varz_chr', chr, '.tsv  ',
                '--make-bed --out  ', num, '__varz_chr', chr
        ), intern=T))
  }

  ## Run merge with PLINK
  merge_list <- cbind(c(paste0(num, '__varz_chr', c(1:22), '.bed')),
                    c(paste0(num, '__varz_chr', c(1:22), '.bim')),
                    c(paste0(num, '__varz_chr', c(1:22), '.fam'))
  )

  write.table(merge_list, file=paste0(num, '__varz_mergelist.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  system(paste0('./plink ',
                '--merge-list  ', num, '__varz_mergelist.tsv  ',
                '--make-bed --out  ', num, '__varz_chrall  '
  ))

  ### Save annot information for REGENIE
  write.table(regenie_annotationfile, file=paste0(num, '__annotationfile_chrall.tsv'), col.names=F, row.names=F, quote=F)
  write.table(regenie_setfile, file=paste0(num, '__setfile_chrall.tsv'), col.names=F, row.names=F, quote=F)
  write.table(c("Mask1 REGENIE"), file=paste0(num, '__maskdef_chrall.tsv'), col.names=F, row.names=F, quote=F)

  ## Run REGENIE for the MAF<0.1% thresholds; keep only the unrel samples
  try(system(paste0('rm  ', num, '__chrall_disease.regenie')))
  try(system(paste0("head ", num, '__regenie_phenofile.tsv')))
  try(system(paste0("head ", num, '__sampleIDs_unrel.tsv')))
  try(system(paste0("head ", num, '__annotationfile_chrall.tsv')))
  try(system(paste0("head ", num, '__annotationfile_chrall.tsv')))
  try(system(paste0("head ", num, '__setfile_chrall.tsv')))
  try(system(paste0("head ", num, '__maskdef_chrall.tsv')))
  try(system(paste0(regenie_path, ' ',
                '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chrall  ',
                '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                '--keep  ', num, '__sampleIDs_unrel.tsv  ',
                '--anno-file  ', num, '__annotationfile_chrall.tsv ',
                '--set-list  ', num, '__setfile_chrall.tsv ',
                '--mask-def  ', num, '__maskdef_chrall.tsv ',
                '--pThresh  0.99  --out ', num, '__chrall'
  ), intern=FALSE))

  ## Process the REGENIE results
  regenie <- fread(paste0(num, '__chrall_disease.regenie'), stringsAsFactors=F, data.table=F)
  regenie <- regenie[which(!grepl("singleton", regenie$ID)), ]
  regenie$ID <- gsub(".Mask1.0.5", "", regenie$ID)
  regenie$firth.n.sample.alt <- round(regenie$A1FREQ * 2 * regenie$N)
  regenie <- regenie[,c("ID", "firth.n.sample.alt", "BETA", "SE", "EXTRA")]
  colnames(regenie)[3:5] <- c("firth.Est", "firth.Est.SE", "firth.failed")
  regenie$firth.cases <- firth.n.cases
  regenie$firth.controls <- firth.n.controls

  regenie_res_tot <- rbind(regenie_res_tot, regenie)

  ######### MAF<0.0001% masks
  #Create set files for regenie; filter PLINK files to needed variants only!
  regenie_setfile <- NULL
  regenie_annotationfile <- NULL
  cat("\t\textracting variant data from PLINK files for MAF<0.001% threshold ...\n")
  for(chr in c(1:22)){
        cat("\t\t\tbusy with chromosome ", chr, "...\n")
        plink_path <- './plink2'
        regenie_path <- './regenie_v3.2.2.gz_x86_64_Linux_mkl'
        plinkfile <- paste0('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c', chr, '_genotype_variant_sample_QCed')
        plinkfile_type <- "pfile"
        max_maf=0.00001
        max_mac=9
        #carz <- NULL
        #phen1 <- phen0

        ## Use the grouping files
        files2_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_lowmem.RData")
        files2_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highmem.RData")
        files2_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highhighmem.RData")
        files2_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_veryhighmem.RData")
        files2_5 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c", chr, "_missense0.5_popmax0.00001_correct.RData")

        f1 <- get(load(files2_1))
        f2 <- get(load(files2_2)) 
        f3 <- get(load(files2_3)) 
        f4 <- get(load(files2_4)) 
        f5 <- get(load(files2_5))
        class(f1$pos) <- class(f2$pos) <- class(f3$pos) <- class(f4$pos) <- class(f5$pos) <- "numeric"
        class(f1$chr) <- class(f2$chr) <- class(f3$chr) <- class(f4$chr) <- class(f5$chr) <- "numeric"
        group <- bind_rows(f1, f2, f3, f4, f5)

        group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
        group$group_id <- sub("_", "__", group$group_id)
        group$group_id <- sub("_gnomAD_POPMAX0.00001", "_POPMAX0.00001", group$group_id)

        ## Filter to the masks we want to test
        group <- group[group$group_id %in% gene_masks, ]
        #group <- group[,c("varid", "alt", "group_id")]
        group$varid <- paste0("chr", group$varid)

        ## Make part of REGENIE grouping files; collect groupings; will use later on
        group$pseudo_annot <- "REGENIE"
        regenie_annotationfile <- rbind(regenie_annotationfile, group[,c("varid", "group_id", "pseudo_annot")])
        set_inter <- NULL
        for(gr in unique(group$group_id)){
                gr_inter <- group[group$group_id==gr, ]
                chromosome <- gr_inter[1,'chr']
                position <- gr_inter[1,'pos']
                collapse <- paste(gr_inter$varid, collapse=",")
                set_inter <- rbind(set_inter, c(gr, chromosome, position, collapse))
        }
        regenie_setfile <- rbind(regenie_setfile, set_inter)
        
        ## Use PLINK2 to filter to a PLINK file that has the needed variants and samples used in the discovery analysis
        write.table(group$varid, file=paste0(num, '__varz_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
        try(system(paste0(plink_path, ' ',
                '--', plinkfile_type, '  ', plinkfile, '  ',
                '--max-maf ', max_maf, '  --max-mac ', max_mac, '  ',
                '--keep  ', num, '__sampleIDs.tsv  ',
                '--extract  ', num, '__varz_chr', chr, '.tsv  ',
                '--make-bed --out  ', num, '__varz_chr', chr
        ), intern=T))
  }

  ## Run merge with PLINK
  merge_list <- cbind(c(paste0(num, '__varz_chr', c(1:22), '.bed')),
                    c(paste0(num, '__varz_chr', c(1:22), '.bim')),
                    c(paste0(num, '__varz_chr', c(1:22), '.fam'))
  )

  write.table(merge_list, file=paste0(num, '__varz_mergelist.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  system(paste0('./plink ',
                '--merge-list  ', num, '__varz_mergelist.tsv  ',
                '--make-bed --out  ', num, '__varz_chrall  '
  ))

  ### Save annot information for REGENIE
  write.table(regenie_annotationfile, file=paste0(num, '__annotationfile_chrall.tsv'), col.names=F, row.names=F, quote=F)
  write.table(regenie_setfile, file=paste0(num, '__setfile_chrall.tsv'), col.names=F, row.names=F, quote=F)
  write.table(c("Mask1 REGENIE"), file=paste0(num, '__maskdef_chrall.tsv'), col.names=F, row.names=F, quote=F)

  ## Run REGENIE for the MAF<0.001% thresholds; keep only the unrel samples
  try(system(paste0('rm  ', num, '__chrall_disease.regenie')))
  try(system(paste0(regenie_path, ' ',
                '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chrall  ',
                '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                '--keep  ', num, '__sampleIDs_unrel.tsv  ',
                '--anno-file  ', num, '__annotationfile_chrall.tsv ',
                '--set-list  ', num, '__setfile_chrall.tsv ',
                '--mask-def  ', num, '__maskdef_chrall.tsv ',
                '--pThresh  0.99  --out ', num, '__chrall'
  ), intern=T))

  ## Process the REGENIE results
  regenie <- fread(paste0(num, '__chrall_disease.regenie'), stringsAsFactors=F, data.table=F)
  regenie <- regenie[which(!grepl("singleton", regenie$ID)), ]
  regenie$ID <- gsub(".Mask1.0.5", "", regenie$ID)
  regenie$firth.n.sample.alt <- round(regenie$A1FREQ * 2 * regenie$N)
  regenie <- regenie[,c("ID", "firth.n.sample.alt", "BETA", "SE", "EXTRA")]
  colnames(regenie)[3:5] <- c("firth.Est", "firth.Est.SE", "firth.failed")
  regenie$firth.cases <- firth.n.cases
  regenie$firth.controls <- firth.n.controls

  regenie_res_tot <- rbind(regenie_res_tot, regenie)

  ######### MAF<1% masks
  #Create set files for regenie; filter PLINK files to needed variants only!
  regenie_setfile <- NULL
  regenie_annotationfile <- NULL
  cat("\t\textracting variant data from PLINK files for MAF<1% threshold ...\n")
  for(chr in c(1:22)){
        cat("\t\t\tbusy with chromosome ", chr, "...\n")
        plink_path <- './plink2'
        regenie_path <- './regenie_v3.2.2.gz_x86_64_Linux_mkl'
        plinkfile <- paste0('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c', chr, '_genotype_variant_sample_QCed')
        plinkfile_type <- "pfile"
        max_maf=0.01
        max_mac=Inf
        #carz <- NULL
        #phen1 <- phen0

        ## Use the grouping files
        files3_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_lowmem.RData")
        files3_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highmem.RData")
        files3_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highhighmem.RData")
        files3_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_veryhighmem.RData")

        f1 <- get(load(files3_1))
        f2 <- get(load(files3_2)) 
        f3 <- get(load(files3_3)) 
        f4 <- get(load(files3_4)) 
        class(f1$pos) <- class(f2$pos) <- class(f3$pos) <- class(f4$pos) <- "numeric"
        class(f1$chr) <- class(f2$chr) <- class(f3$chr) <- class(f4$chr) <- "numeric"
        group <- bind_rows(f1, f2, f3, f4)

        group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
        group$group_id <- sub("_", "__", group$group_id)
        group$group_id <- sub("_gnomAD_POPMAX0.01", "_POPMAX0.01", group$group_id)

        ## Filter to the masks we want to test
        group <- group[group$group_id %in% gene_masks, ]
        #group <- group[,c("varid", "alt", "group_id")]
        group$varid <- paste0("chr", group$varid)

        ## Make part of REGENIE grouping files; collect groupings; will use later on
        group$pseudo_annot <- "REGENIE"
        regenie_annotationfile <- rbind(regenie_annotationfile, group[,c("varid", "group_id", "pseudo_annot")])
        set_inter <- NULL
        for(gr in unique(group$group_id)){
                gr_inter <- group[group$group_id==gr, ]
                chromosome <- gr_inter[1,'chr']
                position <- gr_inter[1,'pos']
                collapse <- paste(gr_inter$varid, collapse=",")
                set_inter <- rbind(set_inter, c(gr, chromosome, position, collapse))
        }
        regenie_setfile <- rbind(regenie_setfile, set_inter)
        
        ## Use PLINK2 to filter to a PLINK file that has the needed variants and samples used in the discovery analysis
        write.table(group$varid, file=paste0(num, '__varz_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
        try(system(paste0(plink_path, ' ',
                '--', plinkfile_type, '  ', plinkfile, '  ',
                '--max-maf ', max_maf, '  --max-mac ', max_mac, '  ',
                '--keep  ', num, '__sampleIDs.tsv  ',
                '--extract  ', num, '__varz_chr', chr, '.tsv  ',
                '--make-bed --out  ', num, '__varz_chr', chr
        ), intern=T))
  }

  ## Run merge with PLINK
  merge_list <- cbind(c(paste0(num, '__varz_chr', c(1:22), '.bed')),
                    c(paste0(num, '__varz_chr', c(1:22), '.bim')),
                    c(paste0(num, '__varz_chr', c(1:22), '.fam'))
  )

  write.table(merge_list, file=paste0(num, '__varz_mergelist.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  system(paste0('./plink ',
                '--merge-list  ', num, '__varz_mergelist.tsv  ',
                '--make-bed --out  ', num, '__varz_chrall  '
  ))

  ### Save annot information for REGENIE
  write.table(regenie_annotationfile, file=paste0(num, '__annotationfile_chrall.tsv'), col.names=F, row.names=F, quote=F)
  write.table(regenie_setfile, file=paste0(num, '__setfile_chrall.tsv'), col.names=F, row.names=F, quote=F)
  write.table(c("Mask1 REGENIE"), file=paste0(num, '__maskdef_chrall.tsv'), col.names=F, row.names=F, quote=F)

  ## Run REGENIE for the MAF<1% thresholds; keep only the unrel samples
  try(system(paste0('rm  ', num, '__chrall_disease.regenie')))
  try(system(paste0(regenie_path, ' ',
                '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chrall  ',
                '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                '--keep  ', num, '__sampleIDs_unrel.tsv  ',
                '--anno-file  ', num, '__annotationfile_chrall.tsv ',
                '--set-list  ', num, '__setfile_chrall.tsv ',
                '--mask-def  ', num, '__maskdef_chrall.tsv ',
                '--pThresh  0.99  --out ', num, '__chrall'
  ), intern=T))

  ## Process the REGENIE results
  regenie <- fread(paste0(num, '__chrall_disease.regenie'), stringsAsFactors=F, data.table=F)
  regenie <- regenie[which(!grepl("singleton", regenie$ID)), ]
  regenie$ID <- gsub(".Mask1.0.5", "", regenie$ID)
  regenie$firth.n.sample.alt <- round(regenie$A1FREQ * 2 * regenie$N)
  regenie <- regenie[,c("ID", "firth.n.sample.alt", "BETA", "SE", "EXTRA")]
  colnames(regenie)[3:5] <- c("firth.Est", "firth.Est.SE", "firth.failed")
  regenie$firth.cases <- firth.n.cases
  regenie$firth.controls <- firth.n.controls

  regenie_res_tot <- rbind(regenie_res_tot, regenie)

  #########################
  # Merge and save results
  ##########################
  cat('\tSaving final results...\n\n')
  rawassoc_res$ID  <- paste0(rawassoc_res$gene, "__", rawassoc_res$mask)
  rawassoc_res  <- merge(rawassoc_res, regenie_res_tot, by="ID", all=T)

  write.table(rawassoc_res, file=paste0('../summary_results_phewas_all_tests_phecode', num, '_with_firths_results.tsv'),
                        col.names=T, row.names=F, quote=F, sep='\t', append=F)

}
#}
