#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
in1=as.character(args[1])
in2=as.character(args[2])
phenotype=as.character(args[3])
genename=as.character(args[4])
score_method=as.character(args[5])
min.mac=as.numeric(args[6])

.libPaths(c("rpackages4_1_3",.libPaths()))

library(data.table)

#git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#git pull --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#source("/medpop/afib/sjurgens/Rscripts/association_source_v2.R")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("UKBB_200KWES_CVD/Cauchy_test.R")

if(score_method == "Score"){
  pval.method <- "Score.pval"
}else if(score_method == "Score.SPA"){
  pval.method <- "SPA.pval"
}

tot <- NULL
load(paste0(in1))
res <- assoc$results
res <- res[which(grepl("freq0.05", rownames(res))), ]
res$mask <- rownames(res)
res$phenotype <- phenotype
res <- res[res$n.alt >= min.mac, ]
tot <- rbind(tot, res)

load(paste0(in2))
res <- assoc$results
res$mask <- rownames(res)
res <- res[which(grepl("freq0.005", rownames(res))), ]
res$phenotype <- phenotype
res <- res[res$n.alt >= min.mac, ]
tot <- rbind(tot, res)

tot <- tidyr::separate(data=tot, col="mask", into=c("gene", "transcript", "variants", "frequency"), sep="_")
write.table(tot, file=paste0(genename, "_", phenotype, '_rawassociation_results_recessive.tsv'), col.names=T, row.names=F, quote=F, sep='\t')

options(stringsAsFactors = F)
totres <- NULL
gres <- NULL
for(transcript in unique(tot$transcript)){
  transres <- NULL
  inter <- tot[tot$phenotype==phenotype & tot$transcript==transcript, ]
  # lof results
  lof <- inter[inter$variants == "hclof", ]
  if(nrow(lof)>1){
    lof_p <- CCT(lof[,pval.method])
  }else if(nrow(lof)==1){
    lof_p <- lof[,pval.method]
  }else{
    lof_p <- NA    
  }
  line <- c(phenotype, transcript, "hclof", lof_p)
  transres <- rbind(transres, line)
  
  # missense results
  lof <- inter[which(grepl("missense", inter$variants) & !grepl("hclof", inter$variants)), ]
  if(nrow(lof)>1){
    lof_p <- CCT(lof[,pval.method])
  }else if(nrow(lof)==1){
    lof_p <- lof[,pval.method]
  }else{
    lof_p <- NA    
  }
  line <- c(phenotype, transcript, "missense", lof_p)
  transres <- rbind(transres, line)
  
  # lof+missense results
  lof <- inter[which(grepl("hclofmissense", inter$variants)), ]
  if(nrow(lof)>1){
    lof_p <- CCT(lof[,pval.method])
  }else if(nrow(lof)==1){
    lof_p <- lof[,pval.method]
  }else{
    lof_p <- NA    
  }
  line <- c(phenotype, transcript, "hclofmissense", lof_p)
  transres <- as.data.frame(rbind(transres, line))
  colnames(transres) <- c("Phenotype", "Transcript", "Variants", "P")
  class(transres$P) <- "numeric"
  trans_p <- CCT(transres$P)
  line <- c(phenotype, transcript, "VARIANTS_COMBINED", trans_p)
  transres <- rbind(transres, line)
  gres <- rbind(gres, transres)
}
class(gres$P) <- "numeric"
gene_p <- CCT(gres[gres$Variants=="VARIANTS_COMBINED", 'P'])
line <- c(phenotype, "TRANSCRIPTS_COMBINED", "VARIANTS_COMBINED", gene_p)
gres <- rbind(gres, line)
totres <- rbind(totres, gres)

write.table(totres, file=paste0(genename, "_", phenotype, '_cauchy_results_reccesive.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
