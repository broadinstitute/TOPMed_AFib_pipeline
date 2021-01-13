# Adapted from:
# 	Author: topmed analysis pipeline, smgogarten
# 	Link: https://github.com/smgogarten/analysis_pipeline/blob/master/R/vcf2gds.R

start_time <- Sys.time()
##### call library
#### if (!requireNamespace("BiocManager", quietly = TRUE))
####     install.packages("BiocManager",repos="https://cloud.r-project.org")
#### BiocManager::install(c('SeqArray'), dependencies=TRUE, clean=TRUE, ask=FALSE, INSTALL_opts='--no-docs --no-demo --byte-compile')
install.packages("digest",repos="https://cloud.r-project.org")

library(SeqArray)

##### call arguments
args <- commandArgs(trailingOnly=T)
bedfile <- args[1]
famfile <- args[2]
bimfile <- args[3]
outfile <- args[4]
ncpus <- as.numeric(args[5])

##### convert vcf to gds
bed.fn<-bedfile
fam.fn<-famfile
bim.fn<-bimfile
out.fn<-paste0(outfile,".gds")
seqBED2GDS(bed.fn, fam.fn, bim.fn, out.fn,parallel=ncpus, verbose=TRUE)

##### complete conversion
end_time <- Sys.time()
difftime<-end_time - start_time
print(paste("Comuptational time of",round(difftime[[1]],digits=2),"secs"))

###### quit R
sessionInfo()
quit("no")
