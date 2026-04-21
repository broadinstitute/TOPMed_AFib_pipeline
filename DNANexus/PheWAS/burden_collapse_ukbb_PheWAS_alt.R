#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
gdsfile=as.character(args[1])
groupfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
outfile=as.character(args[6])
AF.max=as.numeric(args[7])
MAC.max=as.numeric(args[8])
score.method=as.character(args[9])
key.file=as.character(args[10])
key.file.phecodecol=as.character(args[11])
phenum=as.numeric(args[12])


##########################################################
#### RUN burden collaps
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
source("TOPMed_AFib_pipeline/DNANexus/kernell_variance_component_modfied.R")

key <- fread(key.file, stringsAsFactors=F, data.table=F)
#phecodes <- unique(key[,key.file.phecodecol])
#phecode <- phecodes[phenum]
phecode <- key[phenum, key.file.phecodecol]
ancestry <- key[phenum, "Ancestry"]
nullfile <- gsub("PHENUM", phecode, nullfile)
nullfile <- gsub("ANCESTRY", ancestry, nullfile)
outfile <- gsub("PHENUM", phecode, outfile)
outfile <- gsub("ANCESTRY", ancestry, outfile)

perform_burden_collapse(gdsfile=gdsfile,groupfile=groupfile,phenfile=phenfile,ID_col=ID_col,nullfile=nullfile,outfile=outfile, 
				   burden.test=score.method, collapse=FALSE,
				   AF.max=AF.max, MAC.max=MAC.max, use.weights=FALSE
)

sessionInfo()
quit("no")
