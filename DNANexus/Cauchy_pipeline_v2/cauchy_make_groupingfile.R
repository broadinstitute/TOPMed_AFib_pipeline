#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
lof_annotfile=as.character(args[1])
missense_annotfile=as.character(args[2])
genename=as.character(args[3])

.libPaths(c("rpackages4_1_3",.libPaths()))

library(data.table)

#git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#git pull --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#source("/medpop/afib/sjurgens/Rscripts/association_source_v2.R")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("UKBB_200KWES_CVD/Cauchy_test.R")


library(data.table)
load(paste0(lof_annotfile))
#head(group)
if(genename!="ALL"){
    lof_tot <- group[group$group_id==genename, ]
    #dim(lof)
    lof_tot$varid <- paste0(lof_tot$chr, ":", lof_tot$pos, ":", lof_tot$ref, ":", lof_tot$alt)
    #head(lof)
    #table(lof$TranscriptID)
}else{
    lof_tot <- group  
}

load(paste0(missense_annotfile))
if(genename!="ALL"){
    #head(group)
    missense_tot <- group[group$group_id==genename & group$Dtools>=7, ]
    #dim(missense)
    missense_tot$varid <- paste0(missense_tot$chr, ":", missense_tot$pos, ":", missense_tot$ref, ":", missense_tot$alt)
    #table(missense$TranscriptID)
}else{
    missense_tot <- group
}

missense_cutoffs <- c(0.8, 0.6, 0.4, 0.2)
frequency_cutoffs <- c(0.001, 1e-5)

paths <- .libPaths()
genes = unique(c(missense_tot$group_id, lof_tot$group_id))

#install.packages('foreach')
#install.packages('doParallel')
#library(foreach)
#library(doParallel)

cores=detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
cl <- makeCluster(cores[1]-2, outfile='')
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("rpackages4_1_3"))

for(frequency_cutoff in frequency_cutoffs){    
    cat("Busy with freq cutoff", frequency_cutoff, "...\n")
    #group <- foreach(gene=genes, .inorder=FALSE, .combine='rbind') %dopar% {
    for(gene in genes){
        cat("\tBusy with", which(genes==gene), "out of", length(genes), "...\n")
        .libPaths(paths)
        source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
        
        lof <- lof_tot[lof_tot$group_id==gene, ]
        missense <- missense_tot[missense_tot$group_id==gene, ]
        rez_group <- NULL
        for(transcripts in c("ALL", "CANONICAL", c(unique(missense$TranscriptID)))){
            if(transcripts=="CANONICAL"){
                inter_missense1 <- missense[missense$CANONICAL=="YES" & missense$gnomad_POPMAX <= frequency_cutoff, ]
                inter_lof <- lof[lof$CANONICAL=="YES" & lof$gnomad_POPMAX <= frequency_cutoff, ]
            }else if(transcripts=="ALL"){
                inter_lof <- lof[lof$gnomad_POPMAX <= frequency_cutoff, ]
                inter_missense1 <- missense[missense$gnomad_POPMAX <= frequency_cutoff, ]
                inter_missense1 <- inter_missense1[-(which(duplicated(inter_missense1$varid))), ]
            }else{
                inter_missense1 <- missense[missense$CANONICAL!="YES" & missense$TranscriptID==transcripts & missense$gnomad_POPMAX <= frequency_cutoff, ]
                inter_lof <- lof[lof$CANONICAL!="YES" & lof$TranscriptID==transcripts & lof$gnomad_POPMAX <= frequency_cutoff, ]
            }
            for(missense_cutoff in missense_cutoffs){
                inter_missense <- inter_missense1[inter_missense1$Dprop >= missense_cutoff, ]
                inter_lof_missense <- rbind(inter_lof, inter_missense)
                if(nrow(inter_lof_missense)>0){
                    inter_lof_missense$group_id <- paste0(inter_lof_missense$group_id, "_", transcripts, "_hclofmissense", missense_cutoff, "_freq", frequency_cutoff)
                }
                if(nrow(inter_missense)>0){
                    inter_missense$group_id <- paste0(inter_missense$group_id, "_", transcripts, "_missense", missense_cutoff, "_freq", frequency_cutoff)
                }
                rm <- which(duplicated(inter_lof_missense$varid))
                if(length(rm)>0){inter_lof_missense <- inter_lof_missense[-rm, ]}
                rm <- which(duplicated(inter_missense$varid))
                if(length(rm)>0){inter_missense <- inter_missense[-rm, ]}
                rez_group <- rbind(rez_group, inter_lof_missense, inter_missense)
            }
            if(nrow(inter_lof)>0){
                inter_lof$group_id <- paste0(inter_lof$group_id, "_", transcripts, "_hclof_freq", frequency_cutoff)
            }
            rm <- which(duplicated(inter_lof$varid))
            if(length(rm)>0){inter_lof <- inter_lof[-rm, ]}
            rez_group <- rbind(rez_group, inter_lof)
        }
        rez_group
    }
    save(group, file=paste0(genename, '_multiple_groupingfile_v1_freq', frequency_cutoff, '.RData'))
}
stopCluster(cl)


