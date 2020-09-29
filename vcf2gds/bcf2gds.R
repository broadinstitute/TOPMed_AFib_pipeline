# Adapted from:
# 	Author: topmed analysis pipeline, smgogarten
# 	Link: https://github.com/smgogarten/analysis_pipeline/blob/master/R/vcf2gds.R
start_time <- Sys.time()
##### call library
library(SeqArray)

##### call arguments
args <- commandArgs(trailingOnly=T)
bcf <- args[1]
gds_out <- paste(args[2],".gds",sep="")

##### convert vcf to gds
seqBCF2GDS(bcf, gds_out, storage.option="LZMA_RA", verbose=TRUE)

##### complete conversion
end_time <- Sys.time()
difftime<-end_time - start_time
print(paste("Comuptational time of",round(difftime[[1]],digits=2),"secs"))

###### quit R
sessionInfo()
quit("no")
