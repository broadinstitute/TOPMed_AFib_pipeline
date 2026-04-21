#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
keyfile=as.character(args[1])
overvfile=as.character(args[2])
chunk_num=as.numeric(args[3]) # is a phenotype number ranging from 1 to the number of total phenotypes

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
i <- chunk_num

##### Check outfiles exist #####
inter <- inter2 <- NULL
num <- key[i, 'Phecode']
n.cases <- key[i, 'N_cases']
n.controls <- key[i, 'N_controls']
phenoname <- key[i, 'Name']
category <- overv[overv$meaning==phenoname, 'category']
cat('\n\n\nBusy with phenotype', num, 'which is', phenoname, 'and task', i, 'out of 601 tasks...\n\n')

######################
# Result files
######################

cat('\tchecking if results files are present...\n\n')

#MAF<0.1%, mask1 files
system("mkdir mask1")
system("mkdir mask1/lowmem")
system("mkdir mask1/midmem")
system("mkdir mask1/highmem")
### Assumes that chrom 1-22 exist for all mem-levels...
### Might consider using /mnt/ equivalent 
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask1/lowmem/* mask1/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask1/midmem/* mask1/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask1/highmem/* mask1/highmem/"))
mask1files_1 <- paste0("mask1/lowmem/", c(1:22), "_resultsfile.RData")
mask1files_2 <- paste0("mask1/midmem/", c(1:22), "_resultsfile.RData")
mask1files_3 <- paste0("mask1/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask2")
system("mkdir mask2/lowmem")
system("mkdir mask2/midmem")
system("mkdir mask2/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask2/lowmem/* mask2/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask2/midmem/* mask2/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask2/highmem/* mask2/highmem/"))
mask2files_1 <- paste0("mask2/lowmem/", c(1:22), "_resultsfile.RData")
mask2files_2 <- paste0("mask2/midmem/", c(1:22), "_resultsfile.RData")
mask2files_3 <- paste0("mask2/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask3")
system("mkdir mask3/lowmem")
system("mkdir mask3/midmem")
system("mkdir mask3/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask3/lowmem/* mask3/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask3/midmem/* mask3/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask3/highmem/* mask3/highmem/"))
mask3files_1 <- paste0("mask3/lowmem/", c(1:22), "_resultsfile.RData")
mask3files_2 <- paste0("mask3/midmem/", c(1:22), "_resultsfile.RData")
mask3files_3 <- paste0("mask3/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask4")
system("mkdir mask4/lowmem")
system("mkdir mask4/midmem")
system("mkdir mask4/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask4/lowmem/* mask4/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask4/midmem/* mask4/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask4/highmem/* mask4/highmem/"))
mask4files_1 <- paste0("mask4/lowmem/", c(1:22), "_resultsfile.RData")
mask4files_2 <- paste0("mask4/midmem/", c(1:22), "_resultsfile.RData")
mask4files_3 <- paste0("mask4/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask5")
system("mkdir mask5/lowmem")
system("mkdir mask5/midmem")
system("mkdir mask5/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask5/lowmem/* mask5/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask5/midmem/* mask5/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask5/highmem/* mask5/highmem/"))
mask5files_1 <- paste0("mask5/lowmem/", c(1:22), "_resultsfile.RData")
mask5files_2 <- paste0("mask5/midmem/", c(1:22), "_resultsfile.RData")
mask5files_3 <- paste0("mask5/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask6")
system("mkdir mask6/lowmem")
system("mkdir mask6/midmem")
system("mkdir mask6/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask6/lowmem/* mask6/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask6/midmem/* mask6/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask6/highmem/* mask6/highmem/"))
mask6files_1 <- paste0("mask6/lowmem/", c(1:22), "_resultsfile.RData")
mask6files_2 <- paste0("mask6/midmem/", c(1:22), "_resultsfile.RData")
mask6files_3 <- paste0("mask6/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask7")
system("mkdir mask7/lowmem")
system("mkdir mask7/midmem")
system("mkdir mask7/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask7/lowmem/* mask7/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask7/midmem/* mask7/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask7/highmem/* mask7/highmem/"))
mask7files_1 <- paste0("mask7/lowmem/", c(1:22), "_resultsfile.RData")
mask7files_2 <- paste0("mask7/midmem/", c(1:22), "_resultsfile.RData")
mask7files_3 <- paste0("mask7/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask8")
system("mkdir mask8/lowmem")
system("mkdir mask8/midmem")
system("mkdir mask8/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask8/lowmem/* mask8/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask8/midmem/* mask8/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask8/highmem/* mask8/highmem/"))
mask8files_1 <- paste0("mask8/lowmem/", c(1:22), "_resultsfile.RData")
mask8files_2 <- paste0("mask8/midmem/", c(1:22), "_resultsfile.RData")
mask8files_3 <- paste0("mask8/highmem/", c(1:22), "_resultsfile.RData")

system("mkdir mask9")
system("mkdir mask9/lowmem")
system("mkdir mask9/midmem")
system("mkdir mask9/highmem")
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask9/lowmem/* mask9/lowmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask9/midmem/* mask9/midmem/"))
system(paste0("gsutil -m cp gs://path/to/file/", num, "/mask9/highmem/* mask9/highmem/"))
mask9files_1 <- paste0("mask9/lowmem/", c(1:22), "_resultsfile.RData")
mask9files_2 <- paste0("mask9/midmem/", c(1:22), "_resultsfile.RData")
mask9files_3 <- paste0("mask9/highmem/", c(1:22), "_resultsfile.RData")

total_files <- c(mask1files_1, mask1files_2, mask1files_3,
                 mask2files_1, mask2files_2, mask2files_3,
                 mask3files_1, mask3files_2, mask3files_3,
                 mask4files_1, mask4files_2, mask4files_3,
                 mask5files_1, mask5files_2, mask5files_3, 
                 mask6files_1, mask6files_2, mask6files_3, 
                 mask7files_1, mask7files_2, mask7files_3,
                 mask8files_1, mask8files_2, mask8files_3, 
                 mask9files_1, mask9files_2, mask9files_3
                 )


if(!all(file.exists(total_files))){
    print(total_files[which(file.exists(total_files))])
    cat(paste0(" was/were not found.\n"))
    cat("\n\n\n\n\nWARNING: not all required files found!!! Stopping.\n\n\n\n")
}else{
    cat('\tall results files found. Reading in and merging...\n\n')
    cat('\t\tMAF<0.1% hclofnoflag...\n\n')

    rawassoc_res <- NULL
    
    #Mask1 MAF<0.1% LOF; low + med + high
    inter <- summarydata(files=c(mask1files_1, mask1files_2, mask1files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_POPMAX0.001")
    inter$mask <- "hclof_noflag_POPMAX0.001"
    rawassoc_res <- rbind(total, inter)
    
    #MAF<0.1% LOF+missense0.8; low + med + high
    inter <- summarydata(files=c(mask2files_1, mask2files_2, mask2files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_missense0.8_POPMAX0.001")
    inter$mask <- "hclof_noflag_missense0.8_POPMAX0.001"
    rawassoc_res <- rbind(total, inter)

    #MAF<0.1% LOF+missense0.5; low + med + high
    inter <- summarydata(files=c(mask3files_1, mask3files_2, mask3files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_missense0.5_POPMAX0.001")
    inter$mask <- "hclof_noflag_missense0.5_POPMAX0.001"
    rawassoc_res <- rbind(total, inter)

    #MAF<0.001% LOF+missense0.5; low + med + high
    inter <- summarydata(files=c(mask4files_1, mask4files_2, mask4files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_missense0.5_POPMAX0.00001")
    inter$mask <- "hclof_noflag_missense0.5_POPMAX0.00001"
    rawassoc_res <- rbind(total, inter)

    #MAF<0.001% missense0.5; low + med + high
    inter <- summarydata(files=c(mask5files_1, mask5files_2, mask5files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__missense0.5_POPMAX0.00001")
    inter$mask <- "missense0.5_POPMAX0.00001"
    rawassoc_res <- rbind(total, inter)

    #MAF<0.001% missense0.2; low + med + high
    inter <- summarydata(files=c(mask6files_1, mask6files_2, mask6files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__missense0.2_POPMAX0.00001")
    inter$mask <- "missense0.2_POPMAX0.00001"
    rawassoc_res <- rbind(total, inter)
    
    #Mask1 MAF<1% LOF; low + med + high
    inter <- summarydata(files=c(mask7files_1, mask7files_2, mask7files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_POPMAX0.01")
    inter$mask <- "hclof_noflag_POPMAX0.01"
    rawassoc_res <- rbind(total, inter)
    
    #MAF<1% LOF+missense0.8; low + med + high
    inter <- summarydata(files=c(mask8files_1, mask8files_2, mask8files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_missense0.8_POPMAX0.01")
    inter$mask <- "hclof_noflag_missense0.8_POPMAX0.01"
    rawassoc_res <- rbind(total, inter)

    #MAF<1% LOF+missense0.5; low + med + high
    inter <- summarydata(files=c(mask9files_1, mask9files_2, mask9files_3), chrs=c(c(1:22), c(1:22), c(1:22)), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    inter$gene <- paste0(rownames(inter), "__hclof_noflag_missense0.5_POPMAX0.01")
    inter$mask <- "hclof_noflag_missense0.5_POPMAX0.01"
    rawassoc_res <- rbind(total, inter)

    rawassoc_res$cases <- n.cases
    rawassoc_res$controls <- n.controls
                           
    #write.table(rawassoc_res, file=paste0('../summary_results_phewas_all_tests_phecode', num, '.tsv'), col.names=T, row.names=F, quote=F, sep='\t', append=T)

    # Create Cauchy file; run first only for truly rare variant masks; this should work as is!!!!
    cat('\trunning Cauchy combinations for MAF<0.1% and 0.001% masks...\n\n')
    inter <- rawassoc_res
    inter <- inter[inter$n.sample.alt>=20, ]
    inter$gene <- gsub("__.*", "", inter$gene)    
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
    write.table(inter1, file=paste0('summary_results_phewas_cauchy_phecode', num, '.tsv'), col.names=T, row.names=F, quote=F, sep='\t', append=F) # This one is important!

  #############################
  # Build phenotype files; Will depend strongly on how you did the phenotype files!!!
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

      ## THIS ONE IS IMPORTANT: Read in the nullmodel from GENESIS
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
  cat("\t\textracting variant data from PLINK files for MAF<0.1% threshold ...\n")
  regenie <- NULL
  for(chr in c(1:22)){
        regenie_setfile <- NULL
        regenie_annotationfile <- NULL
        cat("\t\t\tbusy with chromosome ", chr, "...\n")
        plink_path <- './plink2'
        regenie_path <- './regenie_v3.2.2.gz_x86_64_Linux_mkl'
        #plinkfile <- paste0('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c', chr, '_genotype_variant_sample_QCed')
        system(paste0("dx download exome-seq:/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c", chr, "_genotype_variant_sample_QCed.*"))
        plinkfile <- paste0('ukb23156_c', chr, '_genotype_variant_sample_QCed')
      
        plinkfile_type <- "pfile"
        max_maf=0.001
        max_mac='100000000'
        #carz <- NULL
        #phen1 <- phen0

        ##Use the grouping files used for the analyses; MAF<0.1% groupings
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask1/lowmem/* mask1/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask1/midmem/* mask1/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask1/highmem/* mask1/highmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask2/lowmem/* mask2/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask2/midmem/* mask2/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask2/highmem/* mask2/highmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask3/lowmem/* mask3/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask3/midmem/* mask3/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask3/highmem/* mask3/highmem/"))
         
        files1_1 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files1_2 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files1_3 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files2_1 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files2_2 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files2_3 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files3_1 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files3_2 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
        files3_3 <- paste0("mask1/lowmem/groupingfile_chr", chr, ".RData")
    
        f1 <- bind_rows(get(load(files1_1)), get(load(files1_2)) , get(load(files1_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_POPMAX0.001")

        f2 <- bind_rows(get(load(files2_1)), get(load(files2_2)) , get(load(files2_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_missense0.8_POPMAX0.001")

        f3 <- bind_rows(get(load(files3_1)), get(load(files3_2)) , get(load(files3_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_missense0.5_POPMAX0.001")
   
        class(f1$pos) <- class(f2$pos) <- class(f3$pos) <- "numeric"
        class(f1$chr) <- class(f2$chr) <- class(f3$chr) <- "numeric"
        group <- bind_rows(f1, f2, f3)

        group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
        rm <- which(duplicated(paste0(group$group_id, "__", group$varid)))
        if(length(rm)>0){group <- group[-rm, ]}
      
        ## Filter to the masks we want to test
        group <- group[group$group_id %in% gene_masks, ]
        if(nrow(group)>0){
            #group <- group[,c("varid", "alt", "group_id")]
            group$varid <- paste0("chr", group$varid)

            ## Make part of REGENIE grouping files; collect groupings; will use later on
            group$pseudo_annot <- "REGENIE"
            #regenie_annotationfile <- rbind(regenie_annotationfile, group[,c("varid", "group_id", "pseudo_annot")])
            regenie_annotationfile <- group[,c("varid", "group_id", "pseudo_annot")]
            set_inter <- NULL
            for(gr in unique(group$group_id)){
                gr_inter <- group[group$group_id==gr, ]
                chromosome <- gr_inter[1,'chr']
                position <- gr_inter[1,'pos']
                collapse <- paste(gr_inter$varid, collapse=",")
                set_inter <- rbind(set_inter, c(gr, chromosome, position, collapse))
            }
            #regenie_setfile <- rbind(regenie_setfile, set_inter)
            regenie_setfile <- set_inter
      
            ## Use PLINK2 to filter to a PLINK file that has the needed variants and samples used in the discovery analysis, so also correct MAF filters can be applied!
            rm <- which(group$varid == "chr9:109080945:T:C") # error variant
            if(length(rm)>0){group <- group[-rm, ]}
            write.table(group$varid, file=paste0(num, '__varz_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            try(system(paste0(plink_path, ' ',
                '--', plinkfile_type, '  ', plinkfile, '  ',
                '--max-maf ', max_maf, '  --max-mac ', max_mac, '  ',
                '--keep  ', num, '__sampleIDs.tsv  ',
                '--extract  ', num, '__varz_chr', chr, '.tsv  ',
                '--make-bed --out  ', num, '__varz_chr', chr
            ), intern=T))
  
            ### Save annot information for REGENIE
            write.table(regenie_annotationfile, file=paste0(num, '__annotationfile_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            write.table(regenie_setfile, file=paste0(num, '__setfile_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            write.table(c("Mask1 REGENIE"), file=paste0(num, '__maskdef_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)

            ## Run REGENIE for the MAF<0.1% thresholds; keep only the unrel samples
            try(system(paste0('rm  ', num, '__chr', chr, '_disease.regenie')))
            try(system(paste0("head ", num, '__regenie_phenofile.tsv')))
            try(system(paste0("head ", num, '__sampleIDs_unrel.tsv')))
            try(system(paste0("head ", num, '__annotationfile_chr', chr, '.tsv')))
            try(system(paste0("head ", num, '__annotationfile_chr', chr, '.tsv')))
            try(system(paste0("head ", num, '__setfile_chr', chr, '.tsv')))
            try(system(paste0("head ", num, '__maskdef_chr', chr, '.tsv')))
            try(system(paste0(regenie_path, ' ',
                '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
                '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                '--keep  ', num, '__sampleIDs_unrel.tsv  ',
                '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
                '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
                '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
                '--pThresh  0.99  --out ', num, '__chr', chr, ' '
          ), intern=FALSE))
          try(system(paste0("rm  ", num, '__annotationfile_chr', chr, '.tsv ')))
          try(system(paste0("rm  ", num, '__setfile_chr', chr, '.tsv')))
          try(system(paste0("rm  ", num, '__maskdef_chr', chr, '.tsv')))
          try(system(paste0("rm  ", num, '__varz_chr', chr, '.*')))

          regenie <- bind_rows(regenie, fread(paste0(num, '__chr', chr, '_disease.regenie'), stringsAsFactors=F, data.table=F))
        }
  }

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
  cat("\t\textracting variant data from PLINK files for MAF<0.001% threshold ...\n")
  regenie <- NULL
  for(chr in c(1:22)){
        regenie_setfile <- NULL
        regenie_annotationfile <- NULL
        cat("\t\t\tbusy with chromosome ", chr, "...\n")
        plink_path <- './plink2'
        regenie_path <- './regenie_v3.2.2.gz_x86_64_Linux_mkl'
        #plinkfile <- paste0('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c', chr, '_genotype_variant_sample_QCed')
        plinkfile <- paste0('ukb23156_c', chr, '_genotype_variant_sample_QCed')

        plinkfile_type <- "pfile"
        max_maf=0.00001
        max_mac=9
        #carz <- NULL
        #phen1 <- phen0

        ## Use the grouping files; MAF<0.0001% files
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask4/lowmem/* mask4/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask4/midmem/* mask4/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask4/highmem/* mask4/highmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask5/lowmem/* mask5/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask5/midmem/* mask5/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask5/highmem/* mask5/highmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask6/lowmem/* mask6/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask6/midmem/* mask6/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask6/highmem/* mask6/highmem/"))
         
        files1_1 <- paste0("mask4/lowmem/groupingfile_chr", chr, ".RData")
        files1_2 <- paste0("mask4/lowmem/groupingfile_chr", chr, ".RData")
        files1_3 <- paste0("mask4/lowmem/groupingfile_chr", chr, ".RData")
        files2_1 <- paste0("mask5/lowmem/groupingfile_chr", chr, ".RData")
        files2_2 <- paste0("mask5/lowmem/groupingfile_chr", chr, ".RData")
        files2_3 <- paste0("mask5/lowmem/groupingfile_chr", chr, ".RData")
        files3_1 <- paste0("mask6/lowmem/groupingfile_chr", chr, ".RData")
        files3_2 <- paste0("mask6/lowmem/groupingfile_chr", chr, ".RData")
        files3_3 <- paste0("mask6/lowmem/groupingfile_chr", chr, ".RData")
    
        f1 <- bind_rows(get(load(files1_1)), get(load(files1_2)) , get(load(files1_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_missense0.5_POPMAX0.00001")

        f2 <- bind_rows(get(load(files2_1)), get(load(files2_2)) , get(load(files2_3))) 
        f1$group_id <- paste0(f1$group_id, "__missense0.5_POPMAX0.00001")

        f3 <- bind_rows(get(load(files3_1)), get(load(files3_2)) , get(load(files3_3))) 
        f1$group_id <- paste0(f1$group_id, "__missense0.2_POPMAX0.00001")
   
        class(f1$pos) <- class(f2$pos) <- class(f3$pos) <- "numeric"
        class(f1$chr) <- class(f2$chr) <- class(f3$chr) <- "numeric"
        group <- bind_rows(f1, f2, f3)
        
        group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
        rm <- which(duplicated(paste0(group$group_id, "__", group$varid)))
        if(length(rm)>0){group <- group[-rm, ]}
      
        ## Filter to the masks we want to test
        group <- group[group$group_id %in% gene_masks, ]
        if(nrow(group)>0){
            #group <- group[,c("varid", "alt", "group_id")]
            group$varid <- paste0("chr", group$varid)

            ## Make part of REGENIE grouping files; collect groupings; will use later on
            group$pseudo_annot <- "REGENIE"
            regenie_annotationfile <- group[,c("varid", "group_id", "pseudo_annot")]
            set_inter <- NULL
            for(gr in unique(group$group_id)){
                gr_inter <- group[group$group_id==gr, ]
                chromosome <- gr_inter[1,'chr']
                position <- gr_inter[1,'pos']
                collapse <- paste(gr_inter$varid, collapse=",")
                set_inter <- rbind(set_inter, c(gr, chromosome, position, collapse))
            }
            regenie_setfile <- set_inter
        
            ## Use PLINK2 to filter to a PLINK file that has the needed variants and samples used in the discovery analysis; also needed to apply correct filters!
            rm <- which(group$varid == "chr9:109080945:T:C") # error variant
            if(length(rm)>0){group <- group[-rm, ]}
            write.table(group$varid, file=paste0(num, '__varz_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            try(system(paste0(plink_path, ' ',
                '--', plinkfile_type, '  ', plinkfile, '  ',
                '--max-maf ', max_maf, '  --max-mac ', max_mac, '  ',
                '--keep  ', num, '__sampleIDs.tsv  ',
                '--extract  ', num, '__varz_chr', chr, '.tsv  ',
                '--make-bed --out  ', num, '__varz_chr', chr
            ), intern=T))
    
            ### Save annot information for REGENIE
            write.table(regenie_annotationfile, file=paste0(num, '__annotationfile_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            write.table(regenie_setfile, file=paste0(num, '__setfile_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            write.table(c("Mask1 REGENIE"), file=paste0(num, '__maskdef_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)

            ## Run REGENIE for the MAF<0.001% thresholds; keep only the unrel samples
            try(system(paste0('rm  ', num, '__chr', chr, '_disease.regenie')))
            try(system(paste0(regenie_path, ' ',
                '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, '  ',
                '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                '--keep  ', num, '__sampleIDs_unrel.tsv  ',
                '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
                '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
                '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
                '--pThresh  0.99  --out ', num, '__chr', chr
            ), intern=FALSE))
            try(system(paste0("rm  ", num, '__annotationfile_chr', chr, '.tsv ')))
            try(system(paste0("rm  ", num, '__setfile_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__maskdef_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__varz_chr', chr, '.*')))

            regenie <- bind_rows(regenie, fread(paste0(num, '__chr', chr, '_disease.regenie'), stringsAsFactors=F, data.table=F))
        }
  }

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
  cat("\t\textracting variant data from PLINK files for MAF<1% threshold ...\n")
  regenie <- NULL
  for(chr in c(1:22)){
        regenie_setfile <- NULL
        regenie_annotationfile <- NULL
        cat("\t\t\tbusy with chromosome ", chr, "...\n")
        plink_path <- './plink2'
        regenie_path <- './regenie_v3.2.2.gz_x86_64_Linux_mkl'
        #plinkfile <- paste0('/mnt/project/exome_450k_plink/merged/genotype_variant_sample_QCed/plink/ukb23156_c', chr, '_genotype_variant_sample_QCed')
        plinkfile <- paste0('ukb23156_c', chr, '_genotype_variant_sample_QCed')
        plinkfile_type <- "pfile"
        max_maf=0.01
        max_mac='100000000'
        #carz <- NULL
        #phen1 <- phen0

        ## Use the grouping files; MAF<1% files
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask7/lowmem/* mask7/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask7/midmem/* mask7/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask7/highmem/* mask7/highmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask8/lowmem/* mask8/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask8/midmem/* mask8/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask8/highmem/* mask8/highmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask9/lowmem/* mask9/lowmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask9/midmem/* mask9/midmem/"))
        system(paste0("gsutil -m cp gs://path/to/groupingfiles/", num, "/mask9/highmem/* mask9/highmem/"))
         
        files1_1 <- paste0("mask7/lowmem/groupingfile_chr", chr, ".RData")
        files1_2 <- paste0("mask7/lowmem/groupingfile_chr", chr, ".RData")
        files1_3 <- paste0("mask7/lowmem/groupingfile_chr", chr, ".RData")
        files2_1 <- paste0("mask8/lowmem/groupingfile_chr", chr, ".RData")
        files2_2 <- paste0("mask8/lowmem/groupingfile_chr", chr, ".RData")
        files2_3 <- paste0("mask8/lowmem/groupingfile_chr", chr, ".RData")
        files3_1 <- paste0("mask9/lowmem/groupingfile_chr", chr, ".RData")
        files3_2 <- paste0("mask9/lowmem/groupingfile_chr", chr, ".RData")
        files3_3 <- paste0("mask9/lowmem/groupingfile_chr", chr, ".RData")
    
        f1 <- bind_rows(get(load(files1_1)), get(load(files1_2)) , get(load(files1_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_POPMAX0.01")

        f2 <- bind_rows(get(load(files2_1)), get(load(files2_2)) , get(load(files2_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_missense0.8_POPMAX0.01")

        f3 <- bind_rows(get(load(files3_1)), get(load(files3_2)) , get(load(files3_3))) 
        f1$group_id <- paste0(f1$group_id, "__hclof_noflag_missense0.5_POPMAX0.01")
   
        class(f1$pos) <- class(f2$pos) <- class(f3$pos) <- "numeric"
        class(f1$chr) <- class(f2$chr) <- class(f3$chr) <- "numeric"
        group <- bind_rows(f1, f2, f3)
        
        group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
        rm <- which(duplicated(paste0(group$group_id, "__", group$varid)))
        if(length(rm)>0){group <- group[-rm, ]}
      
        ## Filter to the masks we want to test
        group <- group[group$group_id %in% gene_masks, ]
        if(nrow(group)>0){
            #group <- group[,c("varid", "alt", "group_id")]
            group$varid <- paste0("chr", group$varid)

            ## Make part of REGENIE grouping files; collect groupings; will use later on
            group$pseudo_annot <- "REGENIE"
            regenie_annotationfile <- group[,c("varid", "group_id", "pseudo_annot")]
            set_inter <- NULL
            for(gr in unique(group$group_id)){
                gr_inter <- group[group$group_id==gr, ]
                chromosome <- gr_inter[1,'chr']
                position <- gr_inter[1,'pos']
                collapse <- paste(gr_inter$varid, collapse=",")
                set_inter <- rbind(set_inter, c(gr, chromosome, position, collapse))
            }
            regenie_setfile <- set_inter
        
            ## Use PLINK2 to filter to a PLINK file that has the needed variants and samples used in the discovery analysis
            rm <- which(group$varid == "chr9:109080945:T:C") # error variant
            if(length(rm)>0){group <- group[-rm, ]}
            write.table(group$varid, file=paste0(num, '__varz_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            try(system(paste0(plink_path, ' ',
                '--', plinkfile_type, '  ', plinkfile, '  ',
                '--max-maf ', max_maf, '  --max-mac ', max_mac, '  ',
                '--keep  ', num, '__sampleIDs.tsv  ',
                '--extract  ', num, '__varz_chr', chr, '.tsv  ',
                '--make-bed --out  ', num, '__varz_chr', chr
            ), intern=T))
    
            ### Save annot information for REGENIE
            write.table(regenie_annotationfile, file=paste0(num, '__annotationfile_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            write.table(regenie_setfile, file=paste0(num, '__setfile_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)
            write.table(c("Mask1 REGENIE"), file=paste0(num, '__maskdef_chr', chr, '.tsv'), col.names=F, row.names=F, quote=F)

            ## Run REGENIE for the MAF<1% thresholds; keep only the unrel samples
            try(system(paste0('rm  ', num, '__chr', chr, '_disease.regenie')))
            try(system(paste0(regenie_path, ' ',
                '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, '  ',
                '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                '--keep  ', num, '__sampleIDs_unrel.tsv  ',
                '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
                '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
                '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
                '--pThresh  0.99  --out ', num, '__chr', chr
            ), intern=FALSE))
            try(system(paste0("rm  ", num, '__annotationfile_chr', chr, '.tsv ')))
            try(system(paste0("rm  ", num, '__setfile_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__maskdef_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__varz_chr', chr, '.*')))
            
            regenie <- bind_rows(regenie, fread(paste0(num, '__chr', chr, '_disease.regenie'), stringsAsFactors=F, data.table=F))
        }
  }
  
  regenie <- regenie[which(!grepl("singleton", regenie$ID)), ]
  regenie$ID <- gsub(".Mask1.0.5", "", regenie$ID)
  regenie$firth.n.sample.alt <- round(regenie$A1FREQ * 2 * regenie$N)
  regenie <- regenie[,c("ID", "firth.n.sample.alt", "BETA", "SE", "EXTRA")]
  colnames(regenie)[3:5] <- c("firth.Est", "firth.Est.SE", "firth.failed")
  regenie$firth.cases <- firth.n.cases
  regenie$firth.controls <- firth.n.controls
  ## can add code to make all results with firth carriers (or cases) <20 --> "MISSING". DON'T REMOVE THE LINE FROM THE DATA. 

  regenie_res_tot <- rbind(regenie_res_tot, regenie)  
  
  #########################
  # Merge and save results
  ##########################
  cat('\tSaving final results...\n\n')
  rawassoc_res$ID  <- paste0(rawassoc_res$gene, "__", rawassoc_res$mask)
  rawassoc_res  <- merge(rawassoc_res, regenie_res_tot, by="ID", all=T)

  ###### Perform check on the output
  if(length(which(rawassoc_res$SPA.converged & is.na(rawassoc_res$firth.Est)))>0){
      cat("WARNING: some results were not re-estimated even though they should have been....\n")
      outfile <- paste0('summary_results_phewas_all_tests_phecode', num, '_with_firths_results_WARNINGS.tsv')
  }else{
      outfile <- paste0('summary_results_phewas_all_tests_phecode', num, '_with_firths_results.tsv')
  }
  write.table(rawassoc_res, file=outfile, col.names=T, row.names=F, quote=F, sep='\t', append=F)
  ###### Move outputs to output directory!
  # All effect size file
  system(paste0("gsutil cp ", outfile, " gs://path/to/collected_results/dir/"))
  # Cauchy file
  system(paste0("gsutil cp summary_results_phewas_cauchy_phecode", num, ".tsv  gs://path/to/collected_results/dir/Cauchy/")
}
#}
