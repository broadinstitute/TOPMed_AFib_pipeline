library(SeqArray)
highqulityvariants<-function(num,bucket){

#### copy gdsfile to notebook space
gdsfile<-paste0("freeze.8.chr",num,".pass_only.phased.01_QC.gds")
copycommand<-paste0("gsutil -m cp ",bucket,"/genotype/freeze8_gds/",gdsfile," ./")
system(copycommand)

#### open gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)

#### call allele counts and allele numbers
AC0 <- seqGetData(gds, "annotation/info/AC")
AC <- AC0$data
AN <- seqGetData(gds, "annotation/info/AN")

#### alternative allele frequency
freq0<-AC/AN

#### total number of individuals
allsamples <- seqGetData(gds, "sample.id")

#### calculate the call rates
callrate0<-(AN/2)/length(allsamples)

#### distribution of call rate
summary(callrate0)

##### length check
cat("Total number of variants with callrate",length(callrate0),"\n")
cat("Total number of variants with freqency",length(freq0),"\n")

##### find gds.id and variant information
allvar <- seqGetData(gds, "variant.id")
allpos <- seqGetData(gds, "position")

alleleinfo<- seqGetData(gds,"allele")
alleleinfo1<-gsub(",",":",alleleinfo)

##### create variant ID
varid<-paste(num,allpos,alleleinfo1,sep=":")
cat("Total number variant with variantid",length(varid),"\n")

##### fiter high quality variants
highqual<-which(freq0 >=0.05 & freq0 <=0.95 & callrate0 >= 0.99)
cat("Total number variant with MAF >=0.05 and callrate 0.99",length(highqual),"\n")

##### find gds and variant id
highqual.gds<-allvar[highqual]
highqual.varid<-varid[highqual]
data0<-data.frame(gdsid=highqual.gds,varid=highqual.varid)

##### export information
write.table(data0,paste0("TOPMed_Freeze8_chr",num,"_Highquality_variants.tsv"),col.names=T,row.names=F,quote=F,sep="\t")

##### copy result to the bucket
com0<-paste0("gsutil cp TOPMed_Freeze8_chr",num,"_Highquality_variants.tsv ",bucket,"/QC/High_Quality_TOPMed/")
system(com0,intern=T)
cat("chromosome",num,"done\n")

seqClose(gds)
}

##### bucket<- Sys.getenv('WORKSPACE_BUCKET')
##### for (num in c(1:20,"X")){
##### highqulityvariants(num,bucket=bucket)
##### }
