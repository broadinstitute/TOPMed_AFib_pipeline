#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
in1=as.character(args[1])
in2=as.character(args[2])
phenotype=as.character(args[3])
genename=as.character(args[4])

library(data.table)

tot <- NULL
load(paste0(outin1))
res <- assoc$results
res <- res[which(grepl("freq0.001", rownames(res))), ]
res$mask <- rownames(res)
res$phenotype <- phenotype
res <- res[res$n.alt >=20, ]
tot <- rbind(tot, res)
load(paste0(outin2))
res <- assoc$results
res$mask <- rownames(res)
rownames(res) <- paste0(rownames(res), " ")
res <- res[which(grepl("freq1e-05", rownames(res))) | which(grepl("freq0 ", rownames(res))), ]
res$phenotype <- phenotype
res <- res[res$n.alt >=20, ]
tot <- rbind(tot, res)

tot <- tidyr::separate(data=tot, col="mask", into=c("gene", "transcript", "variants", "frequency"), sep="_")
write.table(tot, file=paste0(genename, "_", phenotype, '_rawassociation_results.tsv'), col.names=T, row.names=F, quote=F, sep='\t')


### P_cauchy
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    pvals <- pvals[-which(is.na(pvals))]
  }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  
  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

options(stringsAsFactors = F)
totres <- NULL
gres <- NULL
for(transcript in unique(tot$transcript)){
  transres <- NULL
  inter <- tot[tot$phenotype==phenotype & tot$transcript==transcript, ]
  # lof results
  lof <- inter[inter$variants == "hclof", ]
  if(nrow(lof)>1){
    lof_p <- CCT(lof$Score.pval)
  }else if(nrow(lof)==1){
    lof_p <- lof$Score.pval
  }else{
    lof_p <- NA    
  }
  line <- c(phenotype, transcript, "hclof", lof_p)
  transres <- rbind(transres, line)
  
  # missense results
  lof <- inter[which(grepl("missense", inter$variants) & !grepl("hclof", inter$variants)), ]
  if(nrow(lof)>1){
    lof_p <- CCT(lof$Score.pval)
  }else if(nrow(lof)==1){
    lof_p <- lof$Score.pval
  }else{
    lof_p <- NA    
  }
  line <- c(phenotype, transcript, "missense", lof_p)
  transres <- rbind(transres, line)
  
  # lof+missense results
  lof <- inter[which(grepl("hclofmissense", inter$variants)), ]
  if(nrow(lof)>1){
    lof_p <- CCT(lof$Score.pval)
  }else if(nrow(lof)==1){
    lof_p <- lof$Score.pval
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

write.table(totres, file=paste0(genename, "_", phenotype, '_cauchy_results.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
