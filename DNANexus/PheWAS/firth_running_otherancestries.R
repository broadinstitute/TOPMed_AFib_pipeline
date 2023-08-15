#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
keyfile=as.character(args[1])
overvfile=as.character(args[2])
chunk_size=as.numeric(args[3])
chunk_num=as.numeric(args[4])
need_to_run_assocs_file=as.character(args[5])


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

#foreach(i=c(1:nrow(key)), .inorder=FALSE) %dopar% {
#for(i in c(1:nrow(key))){ 
inter <- inter2 <- NULL
num <- key[i, 'Phecode']
n.cases <- key[i, 'N_cases']
n.controls <- key[i, 'N_controls']
phenoname <- key[i, 'Name']
category <- overv[overv$meaning==phenoname, 'category']
cat('\n\n\nBusy with phenotype', num, 'which is', phenoname, 'and task', i, 'out of 535 tasks...\n\n')

## Restrict to needed assocs; should be in format gene__phenotype
need_to_run_assocs <- fread(need_to_run_assocs_file, stringsAsFactors=F, data.table=F)
need_to_run_assocs <- tidyr::separate(need_to_run_assocs, col=assocs, into=c("gene", "phenotype"), sep="__")
need_to_run_assocs <- need_to_run_assocs[need_to_run_assocs$phenotype==phenoname, ]
genes_to_run <- need_to_run_assocs$gene

######################
# Result files
######################

if(length(genes_to_run)<1){
    cat("\n\n\n\n\nWARNING: no assocs remaining to run!!! Stopping.\n\n\n\n")
}else{

  #############################
  # Build phenotype files 
  #############################

  if(reestimate_effects){
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
  write.table(phen0[phen0$POP_tight=="AMR", c("FID", "IID")], file=paste0(num, '__sampleIDs_AMR.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  write.table(phen0[phen0$POP_tight=="AFR", c("FID", "IID")], file=paste0(num, '__sampleIDs_AFR.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  write.table(phen0[phen0$POP_tight=="SAS", c("FID", "IID")], file=paste0(num, '__sampleIDs_SAS.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  write.table(phen0[phen0$POP_tight=="EAS", c("FID", "IID")], file=paste0(num, '__sampleIDs_EAS.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  write.table(phen0[phen0$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), c("FID", "IID")], file=paste0(num, '__sampleIDs_nonEUR.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
  
  cat(paste0("\t\ttotal sample size is: ", nrow(phen0), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1, ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0, ]), "\n\n\n"))

  cat(paste0("\t\tEUR sample size is: ", nrow(phen0[phen0$POP_tight=="EUR", ]), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1 & phen0$POP_tight=="EUR", ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0 & phen0$POP_tight=="EUR", ]), "\n\n\n"))

  cat(paste0("\t\tAMR sample size is: ", nrow(phen0[phen0$POP_tight=="AMR", ]), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1 & phen0$POP_tight=="AMR", ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0 & phen0$POP_tight=="AMR", ]), "\n\n\n"))
  
  cat(paste0("\t\tAFR sample size is: ", nrow(phen0[phen0$POP_tight=="AFR", ]), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1 & phen0$POP_tight=="AFR", ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0 & phen0$POP_tight=="AFR", ]), "\n\n\n"))

  cat(paste0("\t\tEAS sample size is: ", nrow(phen0[phen0$POP_tight=="EAS", ]), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1 & phen0$POP_tight=="EAS", ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0 & phen0$POP_tight=="EAS", ]), "\n\n\n"))

  cat(paste0("\t\tSAS sample size is: ", nrow(phen0[phen0$POP_tight=="SAS", ]), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1 & phen0$POP_tight=="SAS", ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0 & phen0$POP_tight=="SAS", ]), "\n\n\n"))

  cat(paste0("\t\tnonEUR sample size is: ", nrow(phen0[phen0$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), ]), "\n"))
  cat(paste0("\t\tn cases is: ", nrow(phen0[phen0$disease==1 & phen0$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), ]), "\n"))
  cat(paste0("\t\tn controls is: ", nrow(phen0[phen0$disease==0 & phen0$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), ]), "\n\n\n"))

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
  firth.n.cases.EUR <- nrow(phen_unrel[phen_unrel$disease==1 & phen_unrel$POP_tight=="EUR", ])
  firth.n.controls.EUR <- nrow(phen_unrel[phen_unrel$disease==0 & phen_unrel$POP_tight=="EUR", ])    
  firth.n.cases.AMR <- nrow(phen_unrel[phen_unrel$disease==1 & phen_unrel$POP_tight=="AMR", ])
  firth.n.controls.AMR <- nrow(phen_unrel[phen_unrel$disease==0 & phen_unrel$POP_tight=="AMR", ])
  firth.n.cases.AFR <- nrow(phen_unrel[phen_unrel$disease==1 & phen_unrel$POP_tight=="AFR", ])
  firth.n.controls.AFR <- nrow(phen_unrel[phen_unrel$disease==0 & phen_unrel$POP_tight=="AFR", ])
  firth.n.cases.EAS <- nrow(phen_unrel[phen_unrel$disease==1 & phen_unrel$POP_tight=="EAS", ])
  firth.n.controls.EAS <- nrow(phen_unrel[phen_unrel$disease==0 & phen_unrel$POP_tight=="EAS", ])
  firth.n.cases.SAS <- nrow(phen_unrel[phen_unrel$disease==1 & phen_unrel$POP_tight=="SAS", ])
  firth.n.controls.SAS <- nrow(phen_unrel[phen_unrel$disease==0 & phen_unrel$POP_tight=="SAS", ])
  firth.n.cases.nonEUR <- nrow(phen_unrel[phen_unrel$disease==1 & phen_unrel$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), ])
  firth.n.controls.nonEUR <- nrow(phen_unrel[phen_unrel$disease==0 & phen_unrel$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), ])
      
  cat(paste0("\n\n\n\t\tunrelated sample size is: ", nrow(phen_unrel), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls, "\n\n\n"))

  cat(paste0("\t\tEUR unrelated sample size is: ", nrow(phen_unrel[phen_unrel$POP_tight=="EUR", ]), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases.EUR, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls.EUR, "\n\n\n"))

  cat(paste0("\t\tAMR unrelated sample size is: ", nrow(phen_unrel[phen_unrel$POP_tight=="AMR", ]), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases.AMR, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls.AMR, "\n\n\n"))
  
  cat(paste0("\t\tAFR unrelated sample size is: ", nrow(phen_unrel[phen_unrel$POP_tight=="AFR", ]), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases.AFR, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls.AFR, "\n\n\n"))

  cat(paste0("\t\tEAS unrelated sample size is: ", nrow(phen_unrel[phen_unrel$POP_tight=="EAS", ]), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases.EAS, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls.EAS, "\n\n\n"))

  cat(paste0("\t\tSAS unrelated sample size is: ", nrow(phen_unrel[phen_unrel$POP_tight=="SAS", ]), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases.SAS, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls.SAS, "\n\n\n"))

  cat(paste0("\t\tnonEUR unrelated sample size is: ", nrow(phen_unrel[phen_unrel$POP_tight %in% c("AMR", "AFR", "EAS", "SAS"), ]), "\n"))
  cat(paste0("\t\tn cases is: ", firth.n.cases.nonEUR, "\n"))
  cat(paste0("\t\tn controls is: ", firth.n.controls.nonEUR, "\n\n\n"))

      
  ##########################################
  ## Run PLINK -> REGENIE Firth pipeline
  ##########################################

  ###### Identify tests to rerun
  cat('\tRunning PLINK extracting and REGENIE effect size estimation...\n\n')
  ### identify based on whetehr SPA was applied (if SPA was applied, then the nominal P<0.05); also rerun P<0.1 with beta>0
  #assocs <- rawassoc_res[(rawassoc_res$SPA.converged & !is.na(rawassoc_res$SPA.converged)) | (rawassoc_res$SPA.pval<0.1 & rawassoc_res$Est>0), ]
  gene_masks <- NULL
  for(jj in c(1:length(genes_to_run))){
    gene_masks <- c(gene_masks, paste0(genes_to_run[jj], c("__hclofnoflag_POPMAX0.001", 
                                                             "__hclofnoflag_missense0.8_POPMAX0.001",
                                                             "__hclofnoflag_missense0.5_POPMAX0.001",
                                                             "__hclofnoflag_missense0.5_POPMAX0.00001",
                                                             "__missense0.5_POPMAX0.00001",
                                                             "__missense0.2_POPMAX0.00001",
                                                             "__hclofnoflag_POPMAX0.01",
                                                             "__hclofnoflag_missense0.8_POPMAX0.01",
                                                             "__hclofnoflag_missense0.8_POPMAX0.01"))
                     )
  }
  
  regenie_res_tot <- NULL
  regenie_res_tot_nonEUR <- NULL

  ############# FROM HERE: it should be very comparable to the total cohort analyses scripts shard before; differences are highlighted explicitly!  #############
      
  ######### MAF<0.1% masks
  #Create set files for regenie; filter PLINK files to needed variants only!
  cat("\t\textracting variant data from PLINK files for MAF<0.1% threshold ...\n")
  regenie <- NULL
  for(chr in c(1:22)){
      if(chr %in% run_chr_vec){
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

        ##Use the grouping files used for the analyses
        #files1_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_lowmem.RData")
        #files1_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highmem.RData")
        #files1_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highhighmem.RData")
        #files1_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_veryhighmem.RData")
        #files1_5 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_popmax0.001.RData")
        system(paste0("dx download -a -f exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_lowmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highhighmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_veryhighmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_popmax0.001.RData")) 
        files1_1 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_lowmem.RData")
        files1_2 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highmem.RData")
        files1_3 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_highhighmem.RData")
        files1_4 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_popmax0.001_veryhighmem.RData")
        files1_5 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_popmax0.001.RData")
     
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
            
            
            ######### HERE IT IS DIFFERENT! ##########

            
            ## Run REGENIE for the MAF<0.1% thresholds; keep only the unrel samples
            #try(system(paste0(regenie_path, ' ',
            #    '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
            #    '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
            #    '--covarFile  ', num, '__regenie_phenofile.tsv   ',
            #    '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
            #    '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
            #    '--keep  ', num, '__sampleIDs_unrel.tsv  ',
            #    '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
            #    '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
            #    '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
            #    '--pThresh  0.99  --out ', num, '__chr', chr, ' '
            #), intern=FALSE))

            for(ancestry in c("EUR", "AMR", "AFR", "EAS", "SAS", "nonEUR")){
                cat(paste0("\t\trunning Firth for ", ancestry, " ancestry...\n"))
                try(system(paste0(regenie_path, ' ',
                    '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
                    '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                    '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                    '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                    '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                    '--keep  ', num, '__sampleIDs_unrel_', ancestry, '.tsv  ',
                    '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
                    '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
                    '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
                    '--pThresh  0.1  --out ', num, '__chr', chr, '_', ancestry, ' '
                ), intern=FALSE))

                regenie_ancestry_inter <- fread(paste0(num, '__chr', chr, '_', ancestry, '_disease.regenie'), stringsAsFactors=F, data.table=F)
                regenie_ancestry_inter <- regenie_ancestry_inter[which(!grepl("singleton", regenie_ancestry_inter$ID)), ]
                regenie_ancestry_inter$ID <- gsub(".Mask1.0.5", "", regenie_ancestry_inter$ID)
                regenie_ancestry_inter$firth.n.sample.alt <- round(regenie_ancestry_inter$A1FREQ * 2 * regenie_ancestry_inter$N)
                regenie_ancestry_inter <- regenie_ancestry_inter[,c("ID", "firth.n.sample.alt", "P", "BETA", "SE", "EXTRA")]
                colnames(regenie_ancestry_inter)[4:6] <- c("firth.Est", "firth.Est.SE", "firth.failed")
                regenie_ancestry_inter$firth.cases <- get(paste0("firth.cases.", ancestry))
                regenie_ancestry_inter$firth.controls <- get(paste0("firth.controls.", ancestry))
                regenie_ancestry_inter$Ancestry <- ancestry
                
                regenie <- bind_rows(regenie, regenie_ancestry_inter)
            
            }
            try(system(paste0("rm  ", num, '__annotationfile_chr', chr, '.tsv ')))
            try(system(paste0("rm  ", num, '__setfile_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__maskdef_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__varz_chr', chr, '.*')))
        }
     }
  }
  regenie_res_tot <- rbind(regenie_res_tot, regenie)


      
  ######### HERE IT IS THE SAME AGAIN. next difference will be explicitly shown! ##########


  ######### MAF<0.0001% masks
  #Create set files for regenie; filter PLINK files to needed variants only!
  cat("\t\textracting variant data from PLINK files for MAF<0.001% threshold ...\n")
  regenie <- NULL
  for(chr in c(1:22)){
      if(chr %in% run_chr_vec){
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

        ## Use the grouping files
        #files2_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_lowmem.RData")
        #files2_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highmem.RData")
        #files2_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highhighmem.RData")
        #files2_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_veryhighmem.RData")
        #files2_5 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c", chr, "_missense0.5_popmax0.00001_correct.RData")
        system(paste0("dx download -a -f exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_lowmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highhighmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_veryhighmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c", chr, "_missense0.5_popmax0.00001_correct.RData"))
        files2_1 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_lowmem.RData")
        files2_2 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highmem.RData")
        files2_3 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_highhighmem.RData")
        files2_4 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflagmissense0.5_missense0.2_popmax0.00001_veryhighmem.RData")
        files2_5 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_missense0.5_popmax0.00001_correct.RData")

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

           
            ######### HERE IT IS DIFFERENT! ##########
            
            
            ## Run REGENIE for the MAF<0.001% thresholds; keep only the unrel samples
            #try(system(paste0(regenie_path, ' ',
            #    '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
            #    '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
            #    '--covarFile  ', num, '__regenie_phenofile.tsv   ',
            #    '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
            #    '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
            #    '--keep  ', num, '__sampleIDs_unrel.tsv  ',
            #    '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
            #    '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
            #    '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
            #    '--pThresh  0.99  --out ', num, '__chr', chr, ' '
            #), intern=FALSE))

            for(ancestry in c("EUR", "AMR", "AFR", "EAS", "SAS", "nonEUR")){
                cat(paste0("\t\trunning Firth for ", ancestry, " ancestry...\n"))
                try(system(paste0(regenie_path, ' ',
                    '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
                    '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                    '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                    '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                    '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                    '--keep  ', num, '__sampleIDs_unrel_', ancestry, '.tsv  ',
                    '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
                    '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
                    '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
                    '--pThresh  0.1  --out ', num, '__chr', chr, '_', ancestry, ' '
                ), intern=FALSE))

                regenie_ancestry_inter <- fread(paste0(num, '__chr', chr, '_', ancestry, '_disease.regenie'), stringsAsFactors=F, data.table=F)
                regenie_ancestry_inter <- regenie_ancestry_inter[which(!grepl("singleton", regenie_ancestry_inter$ID)), ]
                regenie_ancestry_inter$ID <- gsub(".Mask1.0.5", "", regenie_ancestry_inter$ID)
                regenie_ancestry_inter$firth.n.sample.alt <- round(regenie_ancestry_inter$A1FREQ * 2 * regenie_ancestry_inter$N)
                regenie_ancestry_inter <- regenie_ancestry_inter[,c("ID", "firth.n.sample.alt", "P", "BETA", "SE", "EXTRA")]
                colnames(regenie_ancestry_inter)[4:6] <- c("firth.Est", "firth.Est.SE", "firth.failed")
                regenie_ancestry_inter$firth.cases <- get(paste0("firth.cases.", ancestry))
                regenie_ancestry_inter$firth.controls <- get(paste0("firth.controls.", ancestry))
                regenie_ancestry_inter$Ancestry <- ancestry
                
                regenie <- bind_rows(regenie, regenie_ancestry_inter)
    

            }
            try(system(paste0("rm  ", num, '__annotationfile_chr', chr, '.tsv ')))
            try(system(paste0("rm  ", num, '__setfile_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__maskdef_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__varz_chr', chr, '.*')))
            
        }
     }
  }
  regenie_res_tot <- rbind(regenie_res_tot, regenie)



  ######### HERE IT IS THE SAME AGAIN. next difference will be explicitly shown! ##########


      
  ######### MAF<1% masks
  #Create set files for regenie; filter PLINK files to needed variants only!
  cat("\t\textracting variant data from PLINK files for MAF<1% threshold ...\n")
  regenie <- NULL
  for(chr in c(1:22)){
     if(chr %in% run_chr_vec){
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

        ## Use the grouping files
        #files3_1 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_lowmem.RData")
        #files3_2 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highmem.RData")
        #files3_3 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highhighmem.RData")
        #files3_4 <- paste0("/mnt/project/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_veryhighmem.RData")
        system(paste0("dx download -a -f exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_lowmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highhighmem.RData exome-seq:/sjj/projects/phewas/v1/data/grouping_files/mem_split/ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_veryhighmem.RData"))
        files3_1 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_lowmem.RData")
        files3_2 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highmem.RData")
        files3_3 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highhighmem.RData")
        files3_4 <- paste0("ukbb_phewas_v1_groupingfile_c", chr, "_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_veryhighmem.RData")

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

          ######### HERE IT IS DIFFERENT! ##########

            
            ## Run REGENIE for the MAF<1% thresholds; keep only the unrel samples
            #try(system(paste0(regenie_path, ' ',
            #    '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
            #    '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
            #    '--covarFile  ', num, '__regenie_phenofile.tsv   ',
            #    '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
            #    '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
            #    '--keep  ', num, '__sampleIDs_unrel.tsv  ',
            #    '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
            #    '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
            #    '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
            #    '--pThresh  0.99  --out ', num, '__chr', chr, ' '
            #), intern=FALSE))

            for(ancestry in c("EUR", "AMR", "AFR", "EAS", "SAS", "nonEUR")){
                cat(paste0("\t\trunning Firth for ", ancestry, " ancestry...\n"))
                try(system(paste0(regenie_path, ' ',
                    '--step 2  --bt  --ignore-pred  --bed  ', num, '__varz_chr', chr, ' ',
                    '--firth --approx --firth-se --aaf-bins 0.5  --minMAC  1  ',
                    '--covarFile  ', num, '__regenie_phenofile.tsv   ',
                    '--covarCol ', paste(fixef, collapse=","), '  --catCovarList  ', paste(catCovarList, collapse=","), '  ',
                    '--phenoFile  ', num, '__regenie_phenofile.tsv  --phenoCol disease  ',
                    '--keep  ', num, '__sampleIDs_unrel_', ancestry, '.tsv  ',
                    '--anno-file  ', num, '__annotationfile_chr', chr, '.tsv ',
                    '--set-list  ', num, '__setfile_chr', chr, '.tsv ',
                    '--mask-def  ', num, '__maskdef_chr', chr, '.tsv ',
                    '--pThresh  0.1  --out ', num, '__chr', chr, '_', ancestry, ' '
                ), intern=FALSE))

                regenie_ancestry_inter <- fread(paste0(num, '__chr', chr, '_', ancestry, '_disease.regenie'), stringsAsFactors=F, data.table=F)
                regenie_ancestry_inter <- regenie_ancestry_inter[which(!grepl("singleton", regenie_ancestry_inter$ID)), ]
                regenie_ancestry_inter$ID <- gsub(".Mask1.0.5", "", regenie_ancestry_inter$ID)
                regenie_ancestry_inter$firth.n.sample.alt <- round(regenie_ancestry_inter$A1FREQ * 2 * regenie_ancestry_inter$N)
                regenie_ancestry_inter <- regenie_ancestry_inter[,c("ID", "firth.n.sample.alt", "P", "BETA", "SE", "EXTRA")]
                colnames(regenie_ancestry_inter)[4:6] <- c("firth.Est", "firth.Est.SE", "firth.failed")
                regenie_ancestry_inter$firth.cases <- get(paste0("firth.cases.", ancestry))
                regenie_ancestry_inter$firth.controls <- get(paste0("firth.controls.", ancestry))
                regenie_ancestry_inter$Ancestry <- ancestry
                
                regenie <- bind_rows(regenie, regenie_ancestry_inter)
            
            }
            try(system(paste0("rm  ", num, '__annotationfile_chr', chr, '.tsv ')))
            try(system(paste0("rm  ", num, '__setfile_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__maskdef_chr', chr, '.tsv')))
            try(system(paste0("rm  ", num, '__varz_chr', chr, '.*')))
        }
     }
  }
  regenie_res_tot <- rbind(regenie_res_tot, regenie)


  #########################
  # Merge and save results
  ##########################
  cat('\tSaving final results...\n\n')
  write.table(regenie_res_tot, file=paste0('summary_results_otherancestries_Firth_phecode', num, '.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
      
}
#}
