library(SeqArray)

num=22
gds <- seqOpen(paste0("freeze.8.chr",num,".pass_only.phased.gds"), allow.duplicate=T)

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
length(callrate0)
length(freq0)

allvar <- seqGetData(gds, "variant.id")
allpos <- seqGetData(gds, "position")

alleleinfo<- seqGetData(gds,"allele")
alleleinfo1<-gsub(",",":",alleleinfo)

varid<-paste(num,allpos,alleleinfo1,sep=":")

length(varid)

highqual<-which(freq0 >=0.05 & freq0 <=0.95 & callrate0 >= 0.99)

length(highqual)

highqual.gds<-allvar[highqual]
highqual.varid<-varid[highqual]
data0<-data.frame(gdsid=highqual.gds,varid=highqual.varid)

write.table(data0,paste0("TOPMed_Freeze8_chr",num,"_Highquality_variants.tsv"),col.names=T,row.names=F,quote=F,sep="\t")

bucket<- Sys.getenv('WORKSPACE_BUCKET')
com0<-paste0("gsutil cp TOPMed_Freeze8_chr",num,"_Highquality_variants.tsv ",bucket,"/QC/High_Quality_TOPMed/")
com0
system(com0,intern=T)
