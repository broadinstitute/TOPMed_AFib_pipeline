#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
regenie_outfile=as.character(args[1])
cauchy_outfile=as.character(args[2])
minMAC=as.numeric(args[3])
vcMAXAAF=as.numeric(args[4])
lessthan_vcMAXAAF_remove=as.logical(args[5])
maxMAF_cutoff=as.numeric(args[6])
canonical_only=as.logical(args[7])
keep_singletons=as.logical(args[8])
cMAC_col=as.character(args[9])

library(data.table)
.libPaths(c("rpackages4_1_3",.libPaths()))
source("UKBB_200KWES_CVD/Cauchy_test.R")

message("Input filtering...")

dat <- fread(regenie_outfile, stringsAsFactors = F, data.table=F)
# filter failed tests
dat <- dat[which(is.na(dat$EXTRA) | is.null(dat$EXTRA) | grepl("DF=", dat$EXTRA)), ]
# remove singleton masks
if(is.na(keep_singletons)){
    dat <- dat[which(!grepl("singleton", dat$ALLELE1)), ]
}else if(!keep_singletons){
    dat <- dat[which(!grepl("singleton", dat$ALLELE1)), ]
}
# fix pext coding issues
dat$ID <- gsub("pext0.8", "pext80", dat$ID)
dat$ID <- gsub("pext0.9", "pext90", dat$ID)

# Canonical only filter
if(!is.na(canonical_only)){
    if(canonical_only){
        dat <- dat[which(grepl("canonical", dat$ID) | grepl("CANONICAL", dat$ID)), ]
    }
}

# MAF cutoff filter
if(!is.na(maxMAF_cutoff)){
    if(maxMAF_cutoff<vcMAXAAF){
        message("WARNING: 'maxMAF_cutoff' set lower than 'vcMAXAAF' which will yield strange results. Stopping.")
        stop()
    }else{
        mafs <- dat$ALLELE1
        mafs <- paste0("", gsub(".*\\.", "", mafs))
        mafs[which(grepl("singleton", mafs))] <- "0.0"
        mafs[which(!grepl("e", mafs))] <- paste0("0.", gsub(".*\\.", "", mafs[which(!grepl("e", mafs))]))
        mafs <- as.numeric(mafs)
        dat <- dat[which(mafs<=maxMAF_cutoff), ]
    }
}

if(nrow(dat)==0 | "V2" %in% colnames(dat)){
    cat("\n\n\nNo tests in REGENIE output!! Perhaps no REGENIE tests passing filters.\n\n\n")
    burden <- NULL
    write.table(burden, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
}else{
    dat$TRANSCRIPT_ID <- gsub("\\..*", "", dat$ID)
    print(head(dat$TRANSCRIPT_ID))
    dat$GENE_ID <- gsub("__.*", "", dat$TRANSCRIPT_ID)
    print(head(dat$GENE_ID))
    #head(dat)
    
    #### Merge by mask ####
    burden <- dat[dat$TEST=="ADD", ]
    # Remove results with low MAC; if MAC column provided use that for filtering
    if(!is.na(cMAC_col)){
        burden$N.SAMPLE.ALT  <- burden[, cMAC_col]
    }else{
        burden$N.SAMPLE.ALT <- burden$A1FREQ * burden$N * 2
    }
    burden <- burden[burden$N.SAMPLE.ALT>=minMAC, ]
    
    if(nrow(burden)==0){
        cat("\n\n\nNo tests reaching minor allele count.!!\n\n\n")
        burden <- NULL
        write.table(burden, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
    }else{
        colnamezz <- c("TRANSCRIPT_ID", "GENE_ID", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                            "N", "N.SAMPLE.ALT", "BETA", "SE", "CHISQ", "LOG10P")
        message(colnamezz[which(!colnamezz %in% colnames(burden))])
        burden <- burden[,c("TRANSCRIPT_ID", "GENE_ID", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                            "N", "N.SAMPLE.ALT", "BETA", "SE", "CHISQ", "LOG10P")]
        colnames(burden)[c(11:14)] <- paste0("BURDEN_", colnames(burden)[c(11:14)])
        ACATV <- dat[dat$TEST=="ADD-ACATV", c("ID", "CHISQ", "LOG10P")]
        colnames(ACATV)[c(2:3)] <- paste0("ACATV_", colnames(ACATV)[c(2:3)])
        SKAT <- dat[dat$TEST=="ADD-SKAT", c("ID", "CHISQ", "LOG10P")]
        colnames(SKAT)[c(2:3)] <- paste0("SKAT_", colnames(SKAT)[c(2:3)])
        ACATO <- dat[dat$TEST=="ADD-ACATO", c("ID", "CHISQ", "LOG10P")]
        colnames(ACATO)[c(2:3)] <- paste0("ACATO_", colnames(ACATO)[c(2:3)])
        SBAT <- dat[which(grepl("SBAT", dat$TEST)), c("TRANSCRIPT_ID", "LOG10P")] ### gets added later in pipeline
        burden <- merge(burden, ACATV, by="ID", all.x=T, all.y=F)
        burden <- merge(burden, SKAT, by="ID", all.x=T, all.y=F)
        burden <- merge(burden, ACATO, by="ID", all.x=T, all.y=F)
        #burden <- merge(burden, SBAT, by="ID", all.x=T, all.y=F)
        rm(dat)

        message("Merging by mask groupings...")
        ### Merge by mask groupings, e.g. LOF, missense and LOF+missense
        cauchy <- function(line){
            return(CCT(pvals=line, weights=NULL, log10p=TRUE, ignore0s=FALSE, ignore1s=TRUE))
        }
        
        message("\tLOF masks...")
        lof <- NULL
        uniques <- NULL
        length <- NULL
        try(lof <- cbind(burden[which(!grepl("missense", burden$ALLELE1)), 
                       c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
                     burden[which(!grepl("missense", burden$ALLELE1)),                
                       c(which(grepl("BURDEN_", colnames(burden))),
                         which(grepl("ACATV_", colnames(burden))),
                         which(grepl("SKAT_", colnames(burden))))]
        ))
        try(uniques <- unique(lof$ALLELE1))
        try(length <- length(uniques))
        ## find the mask with largest number - which should encompass all masks/variants! - and run that one first.
        num_rawassocs <- NULL
        if(length>1 | is.null(length)){
            for(unique_num in uniques){
                num_rawassocs <- c(num_rawassocs, nrow(lof[lof$ALLELE1==unique_num, ]))
            }
            max_rawassocs <- which(num_rawassocs==max(num_rawassocs))[1]
            if(max_rawassocs!=1){
                uniques <- c(uniques[max_rawassocs], uniques[-max_rawassocs])
            }
        }
        
        if(length==0 | is.null(length)){
            lof <- NULL
        }else if(length==1){
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[1])){
                lof <- lof[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                lof <- lof[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                              "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(lof)[c(7:ncol(lof))] <- paste0(uniques[1], "_", colnames(lof)[c(7:(ncol(lof)))])
            lof$LOF_cauchy_LOG10P <- apply(X=lof[,which(grepl("LOG10P", colnames(lof)))], MARGIN=1, FUN=cauchy)
            lof <- lof[,-(which(colnames(lof)=="ALLELE1"))]
        }else{
            i<-1
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                lofnew <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                lofnew <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                                          "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(lofnew)[c(7:ncol(lofnew))] <- paste0(uniques[i], "_", colnames(lofnew)[c(7:ncol(lofnew))])
            for(i in c(2:length)){
                if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                    inter <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P")]
                }else{
                    inter <- lof[lof$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
                }
                colnames(inter)[c(2:ncol(inter))] <-  paste0(uniques[i], "_", colnames(inter)[c(2:ncol(inter))])
                lofnew <- merge(lofnew, inter, by="TRANSCRIPT_ID", all=T)
            }
            lof <- lofnew
            lof$LOF_cauchy_LOG10P <- apply(X=lof[,which(grepl("LOG10P", colnames(lof)))], MARGIN=1, FUN=cauchy)
            lof <- lof[,-(which(colnames(lof)=="ALLELE1"))]
        }
            
        message("\tmissense masks...")
        missense <- NULL
        uniques <- NULL
        length <- NULL
        try(missense <- cbind(burden[which(!grepl("LOF", burden$ALLELE1) & !grepl("lof", burden$ALLELE1)), 
                            c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
                     burden[which(!grepl("LOF", burden$ALLELE1) &!grepl("lof", burden$ALLELE1)),                
                            c(which(grepl("BURDEN_", colnames(burden))),
                              which(grepl("ACATV_", colnames(burden))), 
                              which(grepl("SKAT_", colnames(burden))))]
        ))
        try(uniques <- unique(missense$ALLELE1))
        try(length <- length(uniques))
        ## find the mask with largest number - which should encompass all masks/variants! - and run that one first.
        num_rawassocs <- NULL
        if(length>1 | is.null(length)){
            for(unique_num in uniques){
                num_rawassocs <- c(num_rawassocs, nrow(missense[missense$ALLELE1==unique_num, ]))
            }
            max_rawassocs <- which(num_rawassocs==max(num_rawassocs))[1]
            if(max_rawassocs!=1){
                uniques <- c(uniques[max_rawassocs], uniques[-max_rawassocs])
            }
        }
        
        if(length==0 | is.null(length)){
            missense <- NULL
        }else if(length==1){
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[1])){
                missense <- missense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                missense <- missense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                        "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(missense)[c(7:ncol(missense))] <- paste0(uniques[1], "_", colnames(missense)[c(7:(ncol(missense)))])
            missense$missense_cauchy_LOG10P <- apply(X=missense[,which(grepl("LOG10P", colnames(missense)))], MARGIN=1, FUN=cauchy)
            missense <- missense[,-(which(colnames(missense)=="ALLELE1"))]
        }else{
            i<-1
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                missensenew <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                missensenew <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                                                         "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]               
            }
            colnames(missensenew)[c(7:ncol(missensenew))] <- paste0(uniques[i], "_", colnames(missensenew)[c(7:ncol(missensenew))])
            for(i in c(2:length)){
                if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                    inter <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P")]
                }else{
                    inter <- missense[missense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
                }
                colnames(inter)[c(2:ncol(inter))] <-  paste0(uniques[i], "_", colnames(inter)[c(2:ncol(inter))])
                missensenew <- merge(missensenew, inter, by="TRANSCRIPT_ID", all=T)
            }
            missense <- missensenew
            missense$missense_cauchy_LOG10P <- apply(X=missense[,which(grepl("LOG10P", colnames(missense)))], MARGIN=1, FUN=cauchy)
            missense <- missense[,-(which(colnames(missense)=="ALLELE1"))]
        }

        message("\tLOF+missense masks...")
        lofmissense1 <- lofmissense <- NULL
        uniques <- NULL
        length <- NULL
        try(lofmissense <- cbind(burden[which(grepl("LOFmissense", burden$ALLELE1) | grepl("lofmissense", burden$ALLELE1) | grepl("LOF_missense", burden$ALLELE1) | grepl("lof_missense", burden$ALLELE1)), 
                            c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
                     burden[which(grepl("LOFmissense", burden$ALLELE1) | grepl("lofmissense", burden$ALLELE1) | grepl("LOF_missense", burden$ALLELE1) | grepl("lof_missense", burden$ALLELE1)),                
                            c(which(grepl("BURDEN_", colnames(burden))),
                              which(grepl("ACATV_", colnames(burden))),
                              which(grepl("SKAT_", colnames(burden))))]
        ))
        try(uniques <- unique(lofmissense$ALLELE1))
        try(length <- length(uniques))
        ## find the mask with largest number - which should encompass all masks/variants! - and run that one first.
        num_rawassocs <- NULL
        if(length>1 | is.null(length)){
            for(unique_num in uniques){
                num_rawassocs <- c(num_rawassocs, nrow(lofmissense[lofmissense$ALLELE1==unique_num, ]))
            }
            max_rawassocs <- which(num_rawassocs==max(num_rawassocs))[1]
            if(max_rawassocs!=1){
                uniques <- c(uniques[max_rawassocs], uniques[-max_rawassocs])
            }
        }
        if(length==0 | is.null(length)){
            lofmissense1 <- lofmissense <- NULL
        }else if(length==1){
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[1])){
                lofmissense <- lofmissense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                lofmissense <- lofmissense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                              "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(lofmissense)[c(7:ncol(lofmissense))] <- paste0(uniques[1], "_", colnames(lofmissense)[c(7:(ncol(lofmissense)))])
            lofmissense$LOFwithflagmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
            lofmissense <- lofmissense[,-(which(colnames(lofmissense)=="ALLELE1"))]
            lofmissense1 <- lofmissense
        }else{
            i<-1
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                lofmissensenew <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                lofmissensenew <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                                                                 "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(lofmissensenew)[c(7:ncol(lofmissensenew))] <- paste0(uniques[i], "_", colnames(lofmissensenew)[c(7:ncol(lofmissensenew))])
            for(i in c(2:length)){
                if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                    inter <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P")]
                }else{
                    inter <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
                }
                colnames(inter)[c(2:ncol(inter))] <-  paste0(uniques[i], "_", colnames(inter)[c(2:ncol(inter))])
                lofmissensenew <- merge(lofmissensenew, inter, by="TRANSCRIPT_ID", all=T)
            }
            lofmissense <- lofmissensenew
            lofmissense$LOFwithflagmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
            lofmissense <- lofmissense[,-(which(colnames(lofmissense)=="ALLELE1"))]
            lofmissense1 <- lofmissense
        }

        lofmissense <- NULL
        uniques <- NULL
        length <- NULL
        try(lofmissense <- cbind(burden[which(grepl("LOFnoflagmissense", burden$ALLELE1) | grepl("lofnoflagmissense", burden$ALLELE1) | grepl("LOFnoflag_missense", burden$ALLELE1) | grepl("lofnoflag_missense", burden$ALLELE1)), 
                            c("ID", "TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N")],
                     burden[which(grepl("LOFnoflagmissense", burden$ALLELE1) | grepl("lofnoflagmissense", burden$ALLELE1) | grepl("LOFnoflag_missense", burden$ALLELE1) | grepl("lofnoflag_missense", burden$ALLELE1)),                
                            c(which(grepl("BURDEN_", colnames(burden))),
                              which(grepl("ACATV_", colnames(burden))),
                              which(grepl("SKAT_", colnames(burden))))]
        ))
        try(uniques <- unique(lofmissense$ALLELE1))
        try(length <- length(uniques))
        ## find the mask with largest number - which should encompass all masks/variants! - and run that one first.
        num_rawassocs <- NULL
        if(length>1 | is.null(length)){
            for(unique_num in uniques){
                num_rawassocs <- c(num_rawassocs, nrow(lofmissense[lofmissense$ALLELE1==unique_num, ]))
            }
            max_rawassocs <- which(num_rawassocs==max(num_rawassocs))[1]
            if(max_rawassocs!=1){
                uniques <- c(uniques[max_rawassocs], uniques[-max_rawassocs])
            }
        }
        
        if(length==0 | is.null(length)){
            lofmissense <- NULL
        }else if(length==1){
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[1])){
                lofmissense <- lofmissense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                lofmissense <- lofmissense[,c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                              "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(lofmissense)[c(7:ncol(lofmissense))] <- paste0(uniques[1], "_", colnames(lofmissense)[c(7:(ncol(lofmissense)))])
            lofmissense$LOFnoflagmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
        }else{
            i<-1
            if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                lofmissensenew <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", "BURDEN_LOG10P")]
            }else{
                lofmissensenew <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "GENE_ID", "ALLELE1", "CHROM", "GENPOS", "N", 
                                                                                 "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
            }
            colnames(lofmissensenew)[c(7:ncol(lofmissensenew))] <- paste0(uniques[i], "_", colnames(lofmissensenew)[c(7:ncol(lofmissensenew))])
            for(i in c(2:length)){
                if(lessthan_vcMAXAAF_remove & !grepl(paste0(vcMAXAAF), uniques[i])){
                    inter <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P")]
                }else{
                    inter <- lofmissense[lofmissense$ALLELE1==uniques[i], c("TRANSCRIPT_ID", "BURDEN_LOG10P", "ACATV_LOG10P", "SKAT_LOG10P")]
                }
                colnames(inter)[c(2:ncol(inter))] <-  paste0(uniques[i], "_", colnames(inter)[c(2:ncol(inter))])
                lofmissensenew <- merge(lofmissensenew, inter, by="TRANSCRIPT_ID", all=T)
            }
            lofmissense <- lofmissensenew
            lofmissense$LOFnoflagmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
            lofmissense <- lofmissense[,-(which(colnames(lofmissense)=="ALLELE1"))]
            lofmissense <- lofmissense[,c(1, 6:ncol(lofmissense))]
        }
        ### Combine the LOFmissense masks
        if(!is.null(lofmissense1) & !is.null(lofmissense)){
            lofmissense <- merge(lofmissense1, lofmissense, by="TRANSCRIPT_ID", all=T)
            lofmissense$LOFmissense_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("LOG10P", colnames(lofmissense)) & !grepl("_cauchy", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
            lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("LOFwithflagmissense_cauchy_LOG10P", "LOFnoflagmissense_cauchy_LOG10P")))]
        }else if(!is.null(lofmissense1)){
            lofmissense <- lofmissense1
            lofmissense$LOFmissense_cauchy_LOG10P <- lofmissense$LOFwithflagmissense_cauchy_LOG10P 
            lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("LOFwithflagmissense_cauchy_LOG10P", "LOFnoflagmissense_cauchy_LOG10P")))]
        }else if(!is.null(lofmissense)){
            lofmissense$LOFmissense_cauchy_LOG10P <- lofmissense$LOFnoflagmissense_cauchy_LOG10P 
            lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("LOFwithflagmissense_cauchy_LOG10P", "LOFnoflagmissense_cauchy_LOG10P")))]
        }
        
        ############ Merge by transcript ############
        message("Merging by transcript...")
        print(colnames(lofmissense))
        print(colnames(missense[,c(1, 6:ncol(missense))]))
        print(colnames(missense))
        try(lofmissense <- merge(lofmissense, missense[,c(1, 6:ncol(missense))], by="TRANSCRIPT_ID", all=T))
        try(lofmissense <- merge(lofmissense, lof[,c(1, 6:ncol(lof))], by="TRANSCRIPT_ID", all=T))
        ### Add SBAT results
        #SBAT <- dat[which(grepl("SBAT", dat$TEST)), c("TRANSCRIPT_ID", "LOG10P")]
        colnames(SBAT)[2] <- "SBAT_cauchy_LOG10P"
        lofmissense <- merge(lofmissense, SBAT, by="TRANSCRIPT_ID", all=T) 
        #lofmissense <- lofmissense[,c(2, 1, 3:ncol(lofmissense))]
        if(length(which(grepl("_cauchy_LOG10P", colnames(lofmissense))))>1){          
            lofmissense$transcript_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("_cauchy_LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
        }else{
            lofmissense$transcript_cauchy_LOG10P <- lofmissense[,which(grepl("_cauchy_LOG10P", colnames(lofmissense)))]
            cat('size file is', ncol(lofmissense), '...\n')
            head(lofmissense)
        }
        lofmissense$transcript_type <- gsub(".*__", "", lofmissense$TRANSCRIPT_ID)
        lofmissense <- lofmissense[,c(2, 1, 3:5, (ncol(lofmissense)), c(6:(ncol(lofmissense)-1)))]
        rm(lof, missense)
                            
        ############ Merge by gene ############
        message("Merging by gene...")
        uniques <- unique(lofmissense$transcript_type)
        print(head(uniques))
        print(head(lofmissense))
        length <- length(uniques)
        if(length==0 | is.null(length)){
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
            if("TRANSCRIPT_ID" %in% colnames(lofmissense) | "transcript_type" %in% colnames(lofmissense)){lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("TRANSCRIPT_ID", "transcript_type")))]}
            lofmissense$gene_cauchy_LOG10P <- apply(X=lofmissense[,which(grepl("transcript_cauchy_LOG10P", colnames(lofmissense)))], MARGIN=1, FUN=cauchy)
        }else{
            lofmissense <- lofmissense[lofmissense$transcript_type==uniques[1], c(1, 7:(ncol(lofmissense)))]
            if("TRANSCRIPT_ID" %in% colnames(lofmissense) | "transcript_type" %in% colnames(lofmissense)){lofmissense <- lofmissense[,-(which(colnames(lofmissense) %in% c("TRANSCRIPT_ID", "transcript_type")))]}
            colnames(lofmissense)[c(2:(ncol(lofmissense)))] <-  paste0(uniques[1], ":", colnames(lofmissense)[c(2:(ncol(lofmissense)))])
            lofmissense$gene_cauchy_LOG10P <- lofmissense[,which(grepl("transcript_cauchy_LOG10P", colnames(lofmissense)))]
        }

        # Remove columns with no meaningful data (depends on input)
        rm_col <- NULL
        for(i in c(1:(ncol(lofmissense)))){
            if(all(is.na(lofmissense[,i]))){rm_col <- c(rm_col,i)}
        }
        if(!is.null(rm_col)) {lofmissense <- lofmissense[,-rm_col]}

        write.table(lofmissense, file=cauchy_outfile, col.names=T, row.names=F, quote=F, sep='\t')
    }
}
