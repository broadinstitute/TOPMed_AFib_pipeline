########
######## call bucket select
bucket<- Sys.getenv('WORKSPACE_BUCKET')

##### copy 1000G and TOPMed Combined file to folder
outfile<-"freeze.8.chrall.pass_only.phased.pruned.1000G.combined.*"
copycombinedfile<-paste0("gsutil cp ",bucket,"/QC/High_Quality_Combined/",outfile, " ./")
system(copycombinedfile,intern=T)

#######
####### seperate 1000G and TOPMed samples
library(data.table)
famfile<-"freeze.8.chrall.pass_only.phased.pruned.1000G.combined.fam"
famfile<-fread(famfile,header=F,data.table=F)
famfile

#######
####### TOPMed samples
outfile<-"TOPMed_afib_sample.tsv"
topmed_sample<-famfile[grep("^NWD",famfile$V1),c("V1","V2")]
write.table(topmed_sample,outfile,col.names=F,row.names=F,sep="\t",quote=F)

#######
####### TOPMed samples with Plink format
inputfile<-"freeze.8.chrall.pass_only.phased.pruned.1000G.combined"
outputfile<-"freeze.8.chrall.pass_only.phased.pruned"
selecttopmed<-paste0("/tmp/plink/plink --bfile ",inputfile," --keep TOPMed_afib_sample.tsv --keep-allele-order --make-bed --out ",outputfile)
system(selecttopmed,intern=T)

#######
####### 1000G samples with Plink format
inputfile<-"freeze.8.chrall.pass_only.phased.pruned.1000G.combined"
outputfile<-"ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned"
select1000g<-paste0("/tmp/plink/plink --bfile ",inputfile," --remove TOPMed_afib_sample.tsv --keep-allele-order --make-bed --out ",outputfile)
system(select1000g,intern=T)

#######
####### 1000G annotation file
samplefile<-"igsr_samples.tsv"
copysamplefile<-paste0("gsutil cp ",bucket,"/QC/High_Quality_1000G/superpop/",samplefile, " ./")
system(copysamplefile,intern=T)

###### prepare for admixture
###### read 1000G fam file
famfile<-"ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.fam"
popfile<-"ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.pop"
dat <- read.delim(samplefile,as.is=T,header=T)
fam <- read.table(famfile,as.is=T,header=F)
fam$order <- 1:dim(fam)[1]
m1 <- merge(fam,dat,by.x="V1",by.y="Sample.name")
m1 <- m1[order(m1$order),]
write.table(m1$Superpopulation.code,popfile,quote=F,row.names=F,col.names=F,sep="\t")

#######
####### perform admixture
bedfile<-"ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.bed"
stepadmixture<-paste0("/tmp/admixture/admixture ",bedfile," 5 --supervised")
system(stepadmixture,intern=T)

#######
####### copy 1000G learned frequency file for TOPMed input
learnpfile<-"ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.5.P"
topmednput<-"freeze.8.chrall.pass_only.phased.pruned.5.P.in"
copytopmedinput<-paste0("cp ",learnpfile," ",topmednput)
system(copytopmedinput,intern=T)

#######
####### perform admixture again to TOPMed
topmedbed<-"freeze.8.chrall.pass_only.phased.pruned.bed"
stepadmixture<-paste0("/tmp/admixture/admixture -P ",topmedbed," 5")
system(stepadmixture,intern=T)


#######
####### finding ethnic column
q_pop <- read.table("ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.pop",as.is=T,header=F)
q_1kg <- read.table("ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.5.Q",as.is=T,header=F)
q_1kg$V6 <- q_pop$V1

for (i in unique(q_1kg$V6)){
  print(i)
  print(summary(q_1kg[which(q_1kg$V6==i),c(1:5)]))
}

#V1: EUR
#V2: EAS
#V3: AMR
#V4: SAS
#V5: AFR
#######
####### reading admixture result file
dat0<-fread("freeze.8.chrall.pass_only.phased.pruned.5.Q",header=F,sep=" ",data.table=F)
colnames(dat0)<-c("EUR","EAS","AMR","SAS","AFR")

#######
####### 80% is cut-off for ethinic group
dat0$POP<-ifelse(dat0$EUR > 0.80,"EUR",
ifelse(dat0$EAS > 0.80,"EAS",
ifelse(dat0$AMR > 0.80,"AMR",
ifelse(dat0$SAS > 0.80,"SAS",
ifelse(dat0$AFR > 0.80,"AFR",
"UND")))))
table(dat0$POP,useNA="always")

#######
####### read in fam file
ind0<-fread("freeze.8.chrall.pass_only.phased.pruned.fam",header=F,data.table=F)
comb0<-data.frame(NWDID=ind0$V1,POP=dat0$POP)
table(comb0$POP,useNA="always")
write.table(comb0,"freeze.8.chrall.pass_only.phased.pruned.fam.5.pop",col.names=T,row.names=F,quote=F,sep="\t")

#######
####### 1000G admixture result
bucket<- Sys.getenv('WORKSPACE_BUCKET')

##### copy 1000G and TOPMed Combined file to folder
outfile<-"ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.*"
admixtureout<-paste0("gsutil cp ",outfile," ",bucket,"/QC/Admixture/")
system(admixtureout,intern=T)

files<-list.files(pattern="freeze.8.chrall.pass_only.phased.pruned.")
admixtureout2<-files[grep("1000G",files,invert=T)]
paste0("gsutil cp ALL.chrall_GRCh38.genotypes.20170504.cleaned_in_comm.pruned.*
ls freeze.8.chrall.pass_only.phased.pruned.*
