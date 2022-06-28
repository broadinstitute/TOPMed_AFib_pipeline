

library(data.table)
load(paste0(lof_annotfile)

head(group)
lof <- group[group$group_id==gene, ]
#dim(lof)
lof$varid <- paste0(lof$chr, ":", lof$pos, ":", lof$ref, ":", lof$alt)
#head(lof)
#table(lof$TranscriptID)

load(paste0(missense_annotfile))
#head(group)
missense <- group[group$group_id==gene & group$Dtools>=7, ]
#dim(missense)
missense$varid <- paste0(missense$chr, ":", missense$pos, ":", missense$ref, ":", missense$alt)
#table(missense$TranscriptID)

rez_group <- NULL
missense_cutoffs <- c(0.8, 0.6, 0.4, 0.2, 0)
frequency_cutoffs <- c(1e-3, 1e-5, 0)

for(frequency_cutoff in frequency_cutoffs){
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
}

group <- rez_group
save(group, file=paste0(gene, '_multiple_groupingfile_v1.RData'))
