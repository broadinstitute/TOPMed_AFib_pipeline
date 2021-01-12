# Adapted from:
# 	Author: topmed analysis pipeline, smgogarten
# 	Link: https://github.com/smgogarten/analysis_pipeline/blob/master/R/vcf2gds.R
start_time <- Sys.time()
##### call library
library(SeqArray)

##### call arguments
args <- commandArgs(trailingOnly=T)
bedfile <- args[1]
famfile <- args[2]
bimfile <- args[3]
ncpus <- as.numeric(args[4])

##### convert vcf to gds
bed.fn<-bedfile
fam.fn<-famfile
bim.fn<-bimfile
out.fn<-gsub("bed","gds",bedfile)
seqBED2GDS(bed.fn, fam.fn, bim.fn, out.fn,parallel=ncpus, verbose=TRUE)

##### complete conversion
end_time <- Sys.time()
difftime<-end_time - start_time
print(paste("Comuptational time of",round(difftime[[1]],digits=2),"secs"))

###### quit R
sessionInfo()
quit("no")
