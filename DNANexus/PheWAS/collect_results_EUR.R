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

## Read in key file to find the phenotype
key <- fread(keyfile, stringsAsFactors=F, data.table=F)
key <- key[key$Ancestry=="EUR", ]
key <- key[key$N_cases >= 50 & key$N_controls >=50, ]
dim(key)

overv <- fread(overvfile, stringsAsFactors = F, data.table=F)
overv <- overv[overv$included_in_mgb_phewas=="yes", ]

## Find chunk to run
splitz <- split(c(1:nrow(key)), ceiling(seq_along(c(1:nrow(key)))/chunk_size))
key <- key[splitz[[chunk_num]], ]

###### Make omnibus file
line <- t(as.data.frame(c("n.site", "n.alt", "n.sample.alt", "Score", "Score.SE", "Score.Stat", "SPA.pval", "Est", "Est.SE", "PVE", "SPA.converged", "chr", "cMAC", "mpos",
                          "Phenotype", "category", "gene", "mask", "cases", "controls"
)
))
write.table(line, file=paste0('summary_results_phewas_all_tests_largechunk', chunk_num, '.tsv'), col.names=F, row.names=F, quote=F, sep='\t')

single_cauchy_row <- t(as.data.frame(c("gene",
                                       "phenotype", "cases", "controls",
                                       "Score__hclof_noflag_canonical",
                                       "Score.SE__hclof_noflag_canonical",
                                       "Est__hclof_noflag_canonical",
                                       "Est.SE__hclof_noflag_canonical",
                                       "SPA.pval__hclof_noflag_canonical",
                                       "Score__hclof_noflag_missense0.8_canonical",
                                       "Score.SE__hclof_noflag_missense0.8_canonical",
                                       "Est__hclof_noflag_missense0.8_canonical",
                                       "Est.SE__hclof_noflag_missense0.8_canonical",
                                       "SPA.pval__hclof_noflag_missense0.8_canonical",
                                       "Score__missense0.5_canonical",
                                       "Score.SE__missense0.5_canonical",
                                       "Est__missense0.5_canonical",
                                       "Est.SE__missense0.5_canonical",
                                       "SPA.pval__missense0.5_canonical",
                                       "P_cauchy"
)))
colnames(single_cauchy_row) <- single_cauchy_row[1,]
single_cauchy_row <- as.data.frame(single_cauchy_row)
write.table(single_cauchy_row, file=paste0('summary_results_phewas_cauchy_largechunk', chunk_num, '.tsv'), col.names=F, row.names=F, quote=F, sep='\t')
single_cauchy_row[1,c(3:21)] <- 0

##### Check outfiles exist #####
#foreach(i=c(1:nrow(key)), .inorder=FALSE) %dopar% {
for(i in c(1:nrow(key))){
  
  inter <- inter2 <- NULL
  num <- key[i, 'Phecode']
  n.cases <- key[i, 'N_cases']
  n.controls <- key[i, 'N_controls']
  phenoname <- key[i, 'Name']
  category <- overv[overv$meaning==phenoname, 'category']
  cat('\n\n\nBusy with phenotype', num, 'which is', phenoname, 'and task', i, 'out of 535 tasks...\n\n')
  files <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/", num, "_results_chr", c(1:22), "_maf0.001_EUR.RData")
  files2 <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/", num, "_results_chr", c(1:22), "_maf0.00001_EUR.RData")
  if(all(file.exists(files)) & all(file.exists(files2))){
    inter <- summarydata(files=files, chrs=c(1:22), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter <- inter[inter$n.sample.alt>=10, ]
    inter$category <- category
    rownames(inter) <- sub("_", "__", rownames(inter))
    rownames(inter) <- sub("gnomAD_POPMAX0.001", "canonical", rownames(inter))
    inter$gene <- rownames(inter)
    
    inter2 <- summarydata(files=files2, chrs=c(1:22), thre_cMAC=1, add_col=TRUE, add_col_name="phenotype", add_col_value=phenoname)
    inter2 <- inter2[inter2$n.sample.alt>=10, ] #Restrict to results with >=10 carriers here, as we might choose a meta-analysis approach of >=10 carriers in a study and >=20 overall
    inter2$category <- category
    rownames(inter2) <- sub("_", "__", rownames(inter2))
    rownames(inter2) <- sub("gnomAD_POPMAX0.00001", "canonical", rownames(inter2))
    inter2$gene <- rownames(inter2)
    inter <- rbind(inter, inter2)
    inter$mask <- gsub(".*__", "", inter$gene)
    inter$gene <- gsub("__.*", "", inter$gene)
    inter$cases <- n.cases
    inter$controls <- n.controls
    
    write.table(inter, file=paste0('summary_results_phewas_all_tests_largechunk', chunk_num, '.tsv'), col.names=F, row.names=F, quote=F, sep='\t', append=T)
    inter <- inter[inter$n.sample.alt>=20, ] #Restrict to results >=20 carriers here
    inter1 <- inter[inter$mask=="hclof_noflag_canonical", c("phenotype", "cases", "controls", "gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter1)[c(5:9)] <- paste0(colnames(inter1)[c(5:9)], "__hclof_noflag_canonical")
    inter2 <- inter[inter$mask=="hclof_noflag_missense0.8_canonical", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter2)[c(2:6)] <- paste0(colnames(inter2)[c(2:6)], "__hclof_noflag_missense0.8_canonical")
    inter3 <- inter[inter$mask=="missense0.5_canonical", c("gene", "Score", "Score.SE", "Est", "Est.SE", "SPA.pval")]
    colnames(inter3)[c(2:6)] <- paste0(colnames(inter3)[c(2:6)], "__missense0.5_canonical")
    inter1 <- merge(inter1, inter2, by="gene", all=T)
    inter1 <- merge(inter1, inter3, by="gene", all=T)
    inter1$phenotype <- phenoname
    colz <- which(grepl("SPA.pval", colnames(inter1)))
    inter1$P_cauchy  <- apply(inter1[,colz], 1, CCT)
    inter1$cases <- n.cases
    inter1$controls <- n.controls
    write.table(inter1, file=paste0('summary_results_phewas_cauchy_largechunk', chunk_num, '_EUR.tsv'), col.names=F, row.names=F, quote=F, sep='\t', append=T)
  }
}
