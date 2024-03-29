args=(commandArgs(TRUE))
gdsfile=as.character(args[1])
groupfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
outfile=as.character(args[6])

##########################################################
#### RUN hclof_noflag_missense0.9 ExtractKernelStatistics
##########################################################

cat('\nReading in packages for analysis...\n')
source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("TOPMed_AFib_pipeline/DNANexus/kernell_variance_component_modfied.R")
source("TOPMed_AFib_pipeline/DNANexus/ExtractKernelStatistics_error_fixed.R")

kernell_variance_component_v2(gdsfile=gdsfile,groupfile=groupfile,phenfile=phenfile,ID_col=ID_col,nullfile=nullfile,outfile=outfile, test="ExtractKernelStatistics", vc.test="Score", AF.max=0.001, MAC.max=Inf,use.weights=FALSE)

sessionInfo()
quit("no")
