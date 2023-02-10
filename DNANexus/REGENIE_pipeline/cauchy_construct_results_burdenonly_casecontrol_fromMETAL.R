#!/usr/bin/env Rscript

###### Requires a custom METAL meta-analysis from REGENIE burdens, including a cMAC columns ######

#### binary traits
args=(commandArgs(TRUE))
regenie_outfile=as.character(args[1])
cauchy_outfile=as.character(args[2])

library(data.table)
.libPaths(c("rpackages4_1_3",.libPaths()))
system(paste0('git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git'))
source("UKBB_200KWES_CVD/Cauchy_test.R")

dat <- fread(regenie_outfile, stringsAsFactors = F, data.table=F)
## We will remove singleton masks here, as these can introduce false-positives for CMP
rm <- which(grepl("singleton", dat$MarkerName))
if(length(rm)>0){dat <- dat[-rm, ]}
if(nrow(dat)==0 | "V2" %in% colnames(dat)){
    cat("\n\n\nNo tests in REGENIE output!! Perhaps no REGENIE tests passing filters.\n\n\n")
    burden <- NULL
    write.table(burden, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
}else{
    dat$TRANSCRIPT_ID <- gsub("\\..*", "", dat$MarkerName)
    dat$GENE_ID <- gsub("__.*", "", dat$MarkerName)
    #head(dat)
    
    #### Merge by mask: in this case all are burden ####
    burden <- dat
    rm(dat)
    burden <- burden[burden$cMAC>=20, ]
    if(nrow(burden)==0){
        cat("\n\n\nNo tests reaching minor allele count.!!\n\n\n")
        burden <- NULL
        write.table(burden, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
    }else{
        burden <- burden[,c("TRANSCRIPT_ID", "GENE_ID", "CHROM", "GENPOS", "MarkerName", "Allele2", "Allele1", "Freq1",
                            "N", "N_cases", "N_controls", "cMAC", "Effect", "StdErr", "P-value")]
        colnames(burden) <- c("TRANSCRIPT_ID", "GENE_ID", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", 
                              "N", "N_cases", "N_controls", "cMAC", "BETA", "SE", "PVALUE")
        
        ### Merge by mask groupings, e.g. LOF, missense and LOF+missense
        cauchy <- function(line){
            return(CCT(pvals=line, weights=NULL, log10p=FALSE, ignore0s=FALSE, ignore1s=TRUE))
        }
        
        lof <- burden[which(!grepl("missense", burden$ALLELE1)), ]
        uniques <- unique(lof$ALLELE1)
        uniques <- uniques[order(uniques)]
        length <- length(uniques)
        if(length==0 | is.null(length)){
            lof <- NULL
        }else if(length==1){
            lof <- lof[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
            colnames(lof)[c(9:ncol(lof))] <- paste0(uniques[1], "_", colnames(lof)[c(9:(ncol(lof)))])
            lof$LOF_cauchy_PVALUE <- lof[,which(grepl("PVALUE", colnames(lof)))]
        }else{
            i<-1
            lofnew <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
            colnames(lofnew)[c(9:ncol(lofnew))] <- paste0(uniques[i], "_", colnames(lofnew)[c(9:ncol(lofnew))])
            for(i in c(2:length)){
                inter <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "PVALUE")]
                colnames(inter)[c(2)] <-  paste0(uniques[i], "_", colnames(inter)[c(2)])
                lofnew <- merge(lofnew, inter, by="TRANSCRIPT_ID", all=T)
            }
            lof <- lofnew
            lof$LOF_cauchy_PVALUE <- apply(X=lof[,which(grepl("PVALUE", colnames(lof)))], MARGIN=1, FUN=cauchy)
        }
        lof <- lof[,-(which(colnames(lof)=="ALLELE1"))]
        
        missense <- burden[which(!grepl("lof", burden$ALLELE1)), 
                            c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
        uniques <- unique(missense$ALLELE1)
        uniques <- uniques[order(uniques)]
        length <- length(uniques)
        if(length==0 | is.null(length)){
            missense <- NULL
        }else if(length==1){
            missense <- missense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
            colnames(missense)[c(9:ncol(missense))] <- paste0(uniques[1], "_", colnames(missense)[c(9:(ncol(missense)))])
            missense$missense_cauchy_PVALUE <- missense[,which(grepl("PVALUE", colnames(missense)))]
        }else{
            i<-1
            missensenew <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
            colnames(missensenew)[c(9:ncol(missensenew))] <- paste0(uniques[i], "_", colnames(missensenew)[c(9:ncol(missensenew))])
            for(i in c(2:length)){
                inter <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "PVALUE")]
                colnames(inter)[c(2)] <-  paste0(uniques[i], "_", colnames(inter)[c(2)])
                missensenew <- merge(missensenew, inter, by="TRANSCRIPT_ID", all=T)
            }
            missense <- missensenew
            missense$missense_cauchy_PVALUE <- apply(X=missense[,which(grepl("PVALUE", colnames(missense)))], MARGIN=1, FUN=cauchy)
        }
        missense <- missense[,-(which(colnames(missense)=="ALLELE1"))]
        
        lofmissense <- burden[which(grepl("lofmissense", burden$ALLELE1)), 
                            c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
        uniques <- unique(lofmissense$ALLELE1)
        uniques <- uniques[order(uniques)]
        length <- length(uniques)
        if(length==0 | is.null(length)){
            lofmissense <- NULL
        }else if(length==1){
            lofmissense <- lofmissense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
            colnames(lofmissense)[c(9:ncol(lofmissense))] <- paste0(uniques[1], "_", colnames(lofmissense)[c(9:(ncol(lofmissense)))])
            lofmissense$lofmissense_cauchy_PVALUE <- lofmissense[,which(grepl("PVALUE", colnames(lofmissense)))]
        }else{
            i<-1
            lofmissensenew <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "N_cases", "N_controls", "PVALUE")]
            colnames(lofmissensenew)[c(9:ncol(lofmissensenew))] <- paste0(uniques[i], "_", colnames(lofmissensenew)[c(9:ncol(lofmissensenew))])
            for(i in c(2:length)){
                inter <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "PVALUE")]
                colnames(inter)[c(2)] <-  paste0(uniques[i], "_", colnames(inter)[c(2)])
                lofmissensenew <- merge(lofmissensenew, inter, by="TRANSCRIPT_ID", all=T)
            }
            lofmissense <- lofmissensenew
            lofmissense$LOFmissense_cauchy_PVALUE <- apply(X=lofmissense[,which(grepl("PVALUE", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
        }
        lofmissense <- lofmissense[,-(which(colnames(lofmissense)=="ALLELE1"))]
        
        
        ### Merge by transcript ###
        try(lofmissense <- merge(lofmissense, missense[,c(1, 8:ncol(missense))], by="TRANSCRIPT_ID", all=T))
        try(lofmissense <- merge(lofmissense, lof[,c(1, 8:ncol(lof))], by="TRANSCRIPT_ID", all=T))
        if(length(which(grepl("_cauchy_PVALUE", colnames(lofmissense))))>1){          
            lofmissense$transcript_cauchy_PVALUE <- apply(X=lofmissense[,which(grepl("_cauchy_PVALUE", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
        }else{
            lofmissense$transcript_cauchy_PVALUE <- lofmissense[,which(grepl("_cauchy_PVALUE", colnames(lofmissense)))]
            cat('size file is', ncol(lofmissense), '...\n')
            head(lofmissense)
        }
        lofmissense$transcript_type <- gsub(".*__", "", lofmissense$TRANSCRIPT_ID)
        lofmissense$transcript_type <- gsub("0", "", lofmissense$transcript_type)
        lofmissense <- lofmissense[,c(2, 1, 3:7, (ncol(lofmissense)), c(8:(ncol(lofmissense)-1)))]
        rm(lof, missense)
        
        ### Merge by gene ###
        uniques <- unique(lofmissense$transcript_type)
        length <- length(uniques)
        if(length==0 | is.null(length)){
            lofmissense <- NULL
        }else if(length>1){
            i<-1
            lofmissensenew <- lofmissense[lofmissense$transcript_type==uniques[i], ]
            colnames(lofmissensenew)[c(9:ncol(lofmissensenew))] <- paste0(uniques[i], ":", colnames(lofmissensenew)[c(9:ncol(lofmissensenew))])
            for(i in c(2:length)){
                inter <- lofmissense[lofmissense$transcript_type==uniques[i], c(1, 9:(ncol(lofmissense)))]
                colnames(inter)[c(2:(ncol(inter)))] <-  paste0(uniques[i], ":", colnames(inter)[c(2:(ncol(inter)))])
                lofmissensenew <- merge(lofmissensenew, inter, by="GENE_ID", all=T)
            }
            lofmissense <- lofmissensenew
            lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("TRANSCRIPT_ID", "transcript_type")))]
            lofmissense$gene_cauchy_PVALUE <- apply(X=lofmissense[,which(grepl("transcript_cauchy_PVALUE", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
        }else{
            lofmissense <- lofmissense[lofmissense$transcript_type==uniques[1], c(1, 9:(ncol(lofmissense)))]
            lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("TRANSCRIPT_ID", "transcript_type")))]
            colnames(lofmissense)[c(2:(ncol(lofmissense)))] <-  paste0(uniques[1], ":", colnames(lofmissense)[c(2:(ncol(lofmissense)))])
            lofmissense$gene_cauchy_PVALUE <- lofmissense[,which(grepl("transcript_cauchy_PVALUE", colnames(lofmissense)))]
        }
        
        lofmissense <- lofmissense[order(lofmissense$gene_cauchy_PVALUE), ]
        write.table(lofmissense, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
    }
}
