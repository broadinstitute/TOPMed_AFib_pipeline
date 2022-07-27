#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
regenie_outfile=as.character(args[1])
cauchy_outfile=as.character(args[2])

library(data.table)
.libPaths(c("rpackages4_1_3",.libPaths()))
source("UKBB_200KWES_CVD/Cauchy_test.R")

dat <- fread(regenie_outfile, stringsAsFactors = F, data.table=F)

dat$TRANSCRIPT_ID <- gsub("\\..*", "", dat$ID)
dat$GENE_ID <- gsub("__.*", "", dat$TRANSCRIPT_ID)
head(dat)

#### Merge by mask ####
burden <- dat[dat$TEST=="ADD", ]
burden$N.SAMPLE.ALT <- burden$A1FREQ * burden$N * 2
burden <- burden[burden$N.SAMPLE.ALT>=20, ]
burden <- burden[,c("TRANSCRIPT_ID", "GENE_ID", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                    "N", "N.SAMPLE.ALT", "BETA", "SE", "CHISQ", "LOG10P")]
colnames(burden)[c(11:14)] <- paste0("BURDEN_", colnames(burden)[c(11:14)])
ACATV <- dat[dat$TEST=="ADD-ACATV", c("ID", "CHISQ", "LOG10P")]
colnames(ACATV)[c(2:3)] <- paste0("ACATV_", colnames(ACATV)[c(2:3)])
SKAT <- dat[dat$TEST=="SKAT", c("ID", "CHISQ", "LOG10P")]
colnames(SKAT)[c(2:3)] <- paste0("SKAT_", colnames(SKAT)[c(2:3)])
ACATO <- dat[dat$TEST=="ADD-ACATO", c("ID", "CHISQ", "LOG10P")]
colnames(ACATO)[c(2:3)] <- paste0("ACATO_", colnames(ACATO)[c(2:3)])
burden <- merge(burden, ACATV, by="ID", all.x=T, all.y=F)
burden <- merge(burden, SKAT, by="ID", all.x=T, all.y=F)
burden <- merge(burden, ACATO, by="ID", all.x=T, all.y=F)
rm(dat)


### Merge by mask groupings, e.g. LOF, missense and LOF+missense
cauchy <- function(line){
    return(CCT(pvals=line, weights=NULL, log10p=TRUE, ignore0s=FALSE, ignore1s=TRUE))
}
lof <- cbind(burden[which(!grepl("missense", burden$ALLELE1)), 
               c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
             burden[which(!grepl("missense", burden$ALLELE1)),                
               c(which(grepl("BURDEN_", colnames(burden))),
               which(grepl("ACATV_", colnames(burden))), which(grepl("SKAT_", colnames(burden))))]
)
uniques <- unique(lof$ALLELE1)
length <- length(uniques)
if(length==0){
    lof <- NULL
}else if(length==1){
    lof <- lof[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                  "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
    colnames(lof)[c(7:ncol(lof))] <- paste0(uniques[1], "_", colnames(lof)[c(7:(ncol(lof)))])
    lof$LOF_cauchy_LOG10P <- apply(X=lof[,which(grepl("LOG10P", colnames(lof)))], MARGIN=1, FUN=cauchy)
}else{
    i<-1
    lofnew <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                              "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
    colnames(lofnew)[c(7:ncol(lofnew))] <- paste0(uniques[i], "_", colnames(lofnew)[c(7:ncol(lofnew))])
    for(i in c(2:length)){
        inter <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
        colnames(inter)[c(2:4)] <-  paste0(uniques[i], "_", colnames(inter)[c(2:4)])
        lofnew <- merge(lofnew, inter, by="TRANSCRIPT_ID", all=T)
    }
    lof <- lofnew
    lof$LOF_cauchy_LOG10P <- apply(X=lof[,which(grepl("LOG10P", colnames(lof)))], MARGIN=1, FUN=cauchy)
}
lof <- lof[,-(which(colnames(lof)=="ALLELE1"))]

missense <- cbind(burden[which(!grepl("LOF", burden$ALLELE1)), 
                    c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
             burden[which(!grepl("LOF", burden$ALLELE1)),                
                    c(which(grepl("BURDEN_", colnames(burden))),
                      which(grepl("ACATV_", colnames(burden))), which(grepl("SKAT_", colnames(burden))))]
)
uniques <- unique(missense$ALLELE1)
length <- length(uniques)
if(length==0){
    missense <- NULL
}else if(length==1){
    missense <- missense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                  "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
    colnames(missense)[c(7:ncol(missense))] <- paste0(uniques[1], "_", colnames(missense)[c(7:(ncol(missense)))])
    missense$missense_cauchy_LOG10P <- apply(X=missense[,which(grepl("LOG10P", colnames(missense)))], MARGIN=1, FUN=cauchy)
}else{
    i<-1
    missensenew <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                             "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
    colnames(missensenew)[c(7:ncol(missensenew))] <- paste0(uniques[i], "_", colnames(missensenew)[c(7:ncol(missensenew))])
    for(i in c(2:length)){
        inter <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
        colnames(inter)[c(2:4)] <-  paste0(uniques[i], "_", colnames(inter)[c(2:4)])
        missensenew <- merge(missensenew, inter, by="TRANSCRIPT_ID", all=T)
    }
    missense <- missensenew
    missense$missense_cauchy_LOG10P <- apply(X=missense[,which(grepl("LOG10P", colnames(missense)))], MARGIN=1, FUN=cauchy)
}
missense <- missense[,-(which(colnames(missense)=="ALLELE1"))]

lofmissense <- cbind(burden[which(grepl("LOFmissense", burden$ALLELE1)), 
                    c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
             burden[which(grepl("LOFmissense", burden$ALLELE1)),                
                    c(which(grepl("BURDEN_", colnames(burden))),
                      which(grepl("ACATV_", colnames(burden))), which(grepl("SKAT_", colnames(burden))))]
)
uniques <- unique(lofmissense$ALLELE1)
length <- length(uniques)
if(length==0){
    lofmissense <- NULL
}else if(length==1){
    lofmissense <- lofmissense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                  "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
    colnames(lofmissense)[c(7:ncol(lofmissense))] <- paste0(uniques[1], "_", colnames(lofmissense)[c(7:(ncol(lofmissense)))])
    lofmissense$lofmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
}else{
    i<-1
    lofmissensenew <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                             "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
    colnames(lofmissensenew)[c(7:ncol(lofmissensenew))] <- paste0(uniques[i], "_", colnames(lofmissensenew)[c(7:ncol(lofmissensenew))])
    for(i in c(2:length)){
        inter <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
        colnames(inter)[c(2:4)] <-  paste0(uniques[i], "_", colnames(inter)[c(2:4)])
        lofmissensenew <- merge(lofmissensenew, inter, by="TRANSCRIPT_ID", all=T)
    }
    lofmissense <- lofmissensenew
    lofmissense$LOFmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
}
lofmissense <- lofmissense[,-(which(colnames(lofmissense)=="ALLELE1"))]
rm(lof, missense)


### Merge by transcript ###
try(lofmissense <- merge(lofmissense, missense[,c(1, 6:ncol(missense))], by="TRANSCRIPT_ID", all=T))
try(lofmissense <- merge(lofmissense, lof[,c(1, 6:ncol(lof))], by="TRANSCRIPT_ID", all=T))
lofmissense$transcript_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("_cauchy_LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
lofmissense$transcript_type <- gsub(".*__", "", lofmissense$TRANSCRIPT_ID)
lofmissense <- lofmissense[,c(2, 1, 3:5, (ncol(lofmissense)), c(6:(ncol(lofmissense)-1)))]

### Merge by gene ###
uniques <- unique(lofmissense$transcript_type)
length <- length(uniques)
if(length==0){
    lofmissense <- NULL
}else if(length>1){
    i<-1
    lofmissensenew <- lofmissense[lofmissense$transcript_type==uniques[i], ]
    colnames(lofmissensenew)[c(7:ncol(lofmissensenew))] <- paste0(uniques[i], ":", colnames(lofmissensenew)[c(7:ncol(lofmissensenew))])
    for(i in c(2:length)){
        inter <- lofmissense[lofmissense$transcript_type==uniques[i], c(1, 7:(ncol(lofmissense)))]
        colnames(inter)[c(2:(ncol(inter)))] <-  paste0(uniques[i], ":", colnames(inter)[c(2:(ncol(inter)))])
        lofmissensenew <- merge(lofmissensenew, inter, by="GENE_ID", all=T)
    }
    lofmissense <- lofmissensenew
}
lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("TRANSCRIPT_ID", "transcript_type")))]
lofmissense$gene_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("transcript_cauchy_LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)

write.table(lofmissense, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
