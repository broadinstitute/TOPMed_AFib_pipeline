source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("TOPMed_AFib_pipeline/DNANexus/kernell_variance_component_modfied.R")
source("TOPMed_AFib_pipeline/DNANexus/ExtractKernelStatistics_error_fixed.R")
cat('\nReading in packages for analysis...\n')

args=(commandArgs(TRUE))
gdsfile=as.character(args[1])
groupfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
outfile=as.character(args[6])
testtype=as.character(args[7])
afcutoff=as.numeric(args[8])

kernell_variance_component_v3(gdsfile=gdsfile,groupfile=groupfile,phenfile=phenfile,ID_col=ID_col,nullfile=nullfile,outfile=outfile, test="ExtractKernelStatistics", vc.test=testtype, AF.max=afcutoff, MAC.max=Inf,use.weights=FALSE)

sessionInfo()
quit("no")
