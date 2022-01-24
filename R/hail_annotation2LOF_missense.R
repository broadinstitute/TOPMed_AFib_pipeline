library(data.table)
library(jsonlite)

##### read arguments
args=(commandArgs(TRUE))
infile=as.character(args[1])

HCoutfile=paste0(infile,"_hclof.RData")
HCgroupfile=paste0(infile,"_hclof_group.RData")

HCnoflaggroupfile=paste0(infile,"_hclof_noflag_group.RData")

missenseoutfile=paste0(infile,"_missense.RData")
missensegroupfile=paste0(infile,"_missense_group.RData")


#####
#setwd("/medpop/afib/data/PHB/53K_WES/annotation/byChr/")
#infile<-paste0("MGB_53K_IBM_PHB_WES_callset_53K_Jan2022.filtered.annotated.tsv.bgz.chr1.tsv")
#outfile<-paste0("TOPMed_Freeze8_LOF_missense_LOFTEE_HC_chr",num,".RData")
#outgoupfile<-paste0("f8_hclof_chr_",num,"_v2.RData")

dat0<-fread(cmd=paste0("egrep HC ",infile),header=F,data.table=F,sep="\t")

hcextract<-function(x){
test<-fromJSON(x[4])
test2<-test$transcript_consequences
if(!is.null(test2)){
test3<-subset(test2,lof=="HC")
if(nrow(test3)>0){
test4<-data.frame(loc=x[1],allele=x[2],varid=x[3],test3,stringsAsFactors=F)
}else{test4<-NULL}
}else{test4<-NULL}
return(test4)
}

###### extract HC information
res0<-apply(dat0,1,hcextract)
res1<-do.call(rbind,lapply(res0,data.frame))
res2<-data.frame(res1,stringsAsFactors=F,row.names=NULL)
res2$allele<-gsub(",",":",gsub("\\]","",gsub("\\[","",gsub("\"","",res2$allele))))
res2$varid<-paste(res2$loc,res2$allele,sep=":")
save(res2,file=HCoutfile)

######
######
res3<-subset(res2,biotype=="protein_coding")
res3$gvarid<-paste0(res3$varid,res3$gene_id)
res4<-res3
#res4<-res3[!duplicated(res3$gvarid),]

library(tidyr)
res4<-separate(data=res4,col="loc",into=c("chr","pos"),remove=F)
res4$chr<-gsub("chr","",res4$chr)
res4$allele<-gsub("\\[|\\]|\"", "",res4$allele)
res4<-separate(data=res4,col="allele",into=c("ref","alt"),remove=F)

group<-res4[,c("gene_id","chr","pos","ref","alt","transcript_id")]
group<-group[order(res4$chr,res4$pos),]
names(group)[1]<-"group_id"

save(group,file=HCgroupfile)

#########
######### hclof noflags
res0<-get(load(HCoutfile))
res2<-subset(res0,is.na(lof_flags))

######
######
res3<-subset(res2,biotype=="protein_coding")
res3$gvarid<-paste0(res3$varid,res3$gene_id)
res4<-res3
#res4<-res3[!duplicated(res3$gvarid),]

library(tidyr)
res4<-separate(data=res4,col="loc",into=c("chr","pos"),remove=F)
res4$chr<-gsub("chr","",res4$chr)
res4$allele<-gsub("\\[|\\]|\"", "",res4$allele)
res4<-separate(data=res4,col="allele",into=c("ref","alt"),remove=F)

group<-res4[,c("gene_id","chr","pos","ref","alt","transcript_id")]
group<-group[order(res4$chr,res4$pos),]
names(group)[1]<-"group_id"

save(group,file=HCnoflaggroupfile)


#####
##### missense variants
dat0<-fread(cmd=paste0("egrep missense ",infile),header=F,data.table=F,sep="\t")

missenseextract<-function(x){
test<-fromJSON(x[4])
test2<-test$transcript_consequences
if(!is.null(test2)){
test3<-test2[grep("missense",test2$consequence_terms),]
if(nrow(test3)>0){
test4<-data.frame(loc=x[1],allele=x[2],varid=x[3],test3,stringsAsFactors=F)
}else{test4<-NULL}
}else{test4<-NULL}
return(test4)
}

###### extract missense information
res0<-apply(dat0[1:100,],1,missenseextract)
res1<-do.call(rbind,lapply(res0,data.frame))
res2<-data.frame(res1,stringsAsFactors=F,row.names=NULL)
res2$allele<-gsub(",",":",gsub("\\]","",gsub("\\[","",gsub("\"","",res2$allele))))
res2$varid<-paste(res2$loc,res2$allele,sep=":")
save(res2,file=missenseoutfile)

######
######
res3<-subset(res2,biotype=="protein_coding")
res3$gvarid<-paste0(res3$varid,res3$gene_id)
res4<-res3
#res4<-res3[!duplicated(res3$gvarid),]

library(tidyr)
res4<-separate(data=res4,col="loc",into=c("chr","pos"),remove=F)
res4$chr<-gsub("chr","",res4$chr)
res4$allele<-gsub("\\[|\\]|\"", "",res4$allele)
res4<-separate(data=res4,col="allele",into=c("ref","alt"),remove=F)
group<-res4[,c("gene_id","chr","pos","ref","alt","transcript_id")]
group<-group[order(res4$chr,res4$pos),]
names(group)[1]<-"group_id"

save(group,file=missensegroupfile)


sessionInfo()
quit("no")

#####
