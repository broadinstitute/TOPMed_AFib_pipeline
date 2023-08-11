#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
keyfile=as.character(args[1])
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

reestimate_effects <- F

#############################
# Key file and overview files
#############################

## Read in key file to find the phenotype
key <- fread(keyfile, stringsAsFactors=F, data.table=F)
key <- key[key$Ancestry=="ALL", ]
key <- key[key$N_cases >= 50 & key$N_controls >=50, ]
dim(key)

#overv <- fread(overvfile, stringsAsFactors = F, data.table=F)
#overv <- overv[overv$included_in_mgb_phewas=="yes", ]

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
category <- "Ancestry_or_Age_outcome"
cat('\n\n\nBusy with phenotype', num, 'which is', phenoname, 'and task', i, 'out of 5 tasks...\n\n')

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
files <- paste0("dx download -a -f ")
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
    ####if(chunk_num %in% c(1:535)){files <- paste0(files, " ", paste0("exome-seq:/sjj/projects/phewas/v1/results/association/", num, "_results_chr", chr, "_maf0.00001.RData"))}
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
    write.table(rawassoc_res, file=paste0('../summary_results_phewas_all_tests_phecode', num, '.tsv'), col.names=T, row.names=F, quote=F, sep='\t', append=T)

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
    write.table(inter1, file=paste0('../summary_results_phewas_cauchy_phecode', num, '.tsv'), col.names=T, row.names=F, quote=F, sep='\t', append=T)

}
#}
