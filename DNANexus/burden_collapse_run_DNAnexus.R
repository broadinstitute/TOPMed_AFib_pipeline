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

perform_burden_collapse(gdsfile=gdsfile,groupfile=groupfile,phenfile=phenfile,ID_col=ID_col,nullfile=nullfile,outfile=outfile, 
				   burden.test="Score.SPA", collapse=FALSE,
				   AF.max=AF.max, MAC.max=Inf, use.weights=FALSE
)

sessionInfo()
quit("no")
