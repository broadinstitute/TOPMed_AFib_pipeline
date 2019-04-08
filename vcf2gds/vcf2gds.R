# Adapted from:
# 	Author: topmed analysis pipeline, smgogarten
# 	Link: https://github.com/smgogarten/analysis_pipeline/blob/master/R/vcf2gds.R
start_time <- Sys.time()
##### call library
library(SeqArray)

##### call arguments
args <- commandArgs(trailingOnly=T)
vcf <- args[1]
gds_out <- paste(args[2],".gds",sep="")
ncpu <- as.numeric(args[3])

##### Read header2
h <- seqVCF_Header("header2.txt")

##### convert vcf to gds
seqVCF2GDS(vcf, gds_out, header=h, storage.option="LZMA_RA", parallel=ncpu, verbose=TRUE)

##### complete conversion
end_time <- Sys.time()
difftime<-end_time - start_time
print(paste("Comuptational time of",round(difftime[[1]],digits=2),"secs"))

###### quit R
sessionInfo()
quit("no")
