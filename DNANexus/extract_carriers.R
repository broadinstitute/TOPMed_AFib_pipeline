#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
groupfile=as.character(args[1])
pfile=as.character(args[2])
plink_path=as.character(args[3])
collapse=as.logical(args[4])
canonical=as.logical(args[5])
outfile=as.character(args[6])
max_maf=as.numeric(args[7])

##########################################################
#### RUN extraction
##########################################################

cat('\nReading in packages for analysis...\n')
### library(GENESIS)
### library(data.table)
### library(SeqArray)
### library(SeqVarTools)

.libPaths(c("rpackages4_1_3",.libPaths()))

#git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#git pull --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#source("/medpop/afib/sjurgens/Rscripts/association_source_v2.R")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")

group <- get(load(groupfile))
group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)

if(canonical){
    group <- group[which(grepl("CANONICAL", group$group_id)), ]
}
group <- group[,c("varid", "alt", "group_id")]
group$varid <- paste0("chr", group$varid)


final <- NULL
num <- 1
length_nums <- length(unique(group$group_id))
for(grouping in unique(group$group_id)){
    cat('busy with number', num, 'out of', length_nums, '...\n')
    num <- num+1
    write.table(group[group$group_id==grouping, 'varid'], file=paste0('varz_', grouping, '_freq', max_maf, '.tsv'), col.names=F, row.names=F, quote=F)
    write.table(group[group$group_id==grouping,c("varid", "alt")], file=paste0('export-allele_', grouping, '_freq', max_maf, '.tsv'), col.names=F, row.names=F, quote=F)
    system(paste0(plink_path, ' ',
                  '--pfile  ', pfile, '  ',
                  '--extract varz_', grouping, '_freq', max_maf, '.tsv  ',
                  '--make-bed --out bfile_', grouping, '_freq', max_maf
    ))
    system(paste0(plink_path, ' ',
                  ' --bfile  bfile_', grouping, '_freq', max_maf,
                  ' --max-maf ', max_maf, ' ',
                  ' --export A --export-allele export-allele_', grouping, "_freq", max_maf, '.tsv',
                  ' --out text_', grouping, '_freq', max_maf
    ))
    library(data.table)
    library(dplyr)
    raw <- fread(paste0('text_', grouping, '_freq', max_maf, '.raw'), stringsAsFactors=F, data.table=F)
    raw <- raw %>% replace(is.na(.), 0)
    
    if(ncol(raw)==6){
        raw[,paste0(grouping)] <- 0
    }else if(ncol(raw)==7){
        raw[,paste0(grouping)] <- raw[,7]
    }else{
        raw[,paste0(grouping)] <- rowSums(raw[,c(7:(ncol(raw)))])
    }
    if(collapse){
        raw[which(raw[,paste0(grouping)]>1), paste0(grouping)] <-1
    }
    colnames(raw)[(ncol(raw))] <- paste0(grouping, "__freq", max_maf) 
    if(is.null(final)){
        raw <- raw[,c(1:6, (ncol(raw)))]
        final <- raw
    }else{
        raw <- raw[,c(1, (ncol(raw)))]
        final <- merge(final, raw, by="FID", all=T)
    }
    system(paste0('rm bfile_', grouping, '_freq', max_maf, '.*'))
    system(paste0('rm text_', grouping, '_freq', max_maf, '.*'))
    
}

write.table(final, file=outfile, col.names=T, row.names=F, quote=F, sep='\t')
