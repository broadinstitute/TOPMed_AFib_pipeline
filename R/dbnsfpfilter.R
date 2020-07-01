#########

######### weighted by tissue expression
######### deleteriousness
##### examples
library(data.table)

dbnsfpfilter<-function(filename,normvar,scorevar=NULL,phredvar=NULL){

name1<-fread(cmd=paste0("zcat ",filename,"| head -n 2"),header=T,data.table=F,sep="\t")
name2<-names(name1)
colnumbers<-which(name2 %in% c(normvar,scorevar,phredvar))

dat1<-fread(cmd=paste0("zcat ",filename),header=T,data.table=F,sep="\t",select=colnumbers)

scorenewvars<-gsub("rankscore","pred",scorevar)
if(length(scorenewvars)>0){

for (score in 1:length(scorenewvars)){
score1<-scorevar[score]
newvar<-scorenewvars[score]
dat1[,newvar]<-ifelse(is.na(dat1[,score1]),NA,ifelse(dat1[,score1]>0.90,"D","B"))
}
}

phrednewvars<-gsub("phred","pred",phredvar)
if(length(phrednewvars)>0){

for (phred in 1:length(phrednewvars)){
phred1<-phredvar[phred]
newvar<-phrednewvars[phred]
dat1[,newvar]<-ifelse(is.na(dat1[,phred1]),NA,ifelse(dat1[,phred1]>0.90,"D","B"))
}
}
dat2<-dat1[,c(normvar,scorenewvars,phrednewvars)]

return(dat2)
}

###### Example 1)
###### filename<-"/medpop/afib/data/NHLBI_WGS/freeze.8/Annotation/WGSA_parsed/v1/chr21_dbNSFP_scatter_piece.gathered.sorted.tsv.gz"
###### normvar<-c("chr","pos","ref","alt","Ensembl_geneid","Ensembl_transcriptid","MetaLR_pred","MetaSVM_pred","MutationAssessor_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred",
###### "PROVEAN_pred","M_CAP_pred","FATHMM_pred","SIFT4G_pred","SIFT_pred","PrimateAI_pred","DEOGEN2_pred")
###### scorevar<-c("MutPred_rankscore","REVEL_rankscore","VEST4_rankscore")
###### res0<-dbsfpfilter(filename=filename,normvar=normvar,scorevar=scorevar)


###### Example 2)
###### num=21
###### filename<-paste0("/medpop/afib/data/NHLBI_WGS/freeze.8/Annotation/WGSA_parsed/v1/chr",num,"_VEP_scatter_piece.tsv.gathered.sorted.tsv.gz")
###### normvar<-c("chr","pos","ref","alt","VEP_ensembl_Gene_ID","VEP_ensembl_Gene_Name","VEP_ensembl_Transcript_ID","VEP_ensembl_Consequence","VEP_ensembl_LoF")
###### scorevar<-c("DANN_rank_score","fathmm_MKL_coding_rankscore","fathmm_XF_rankscore","GenoCanyon_rankscore")
###### phredvar<-c("CADD_phred","Eigen_phred","Eigen_PC_phred")
###### outfile<-paste0("/medpop/afib/data/NHLBI_WGS/freeze.8/Annotation/WGSA_parsed/v1/chr",num,"_VEP_scatter_piece.tsv.gathered.sorted.pred.tsv")





####### dbnsfp dataset
####### filename<-paste0("/medpop/afib/data/NHLBI_WGS/freeze.8/Annotation/WGSA_parsed/v1/chr",num,"_dbNSFP_scatter_piece.gathered.sorted.tsv.gz")
####### normvar<-c("chr","pos","ref","alt","Ensembl_geneid","Ensembl_transcriptid","MetaLR_pred","MetaSVM_pred","MutationAssessor_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred",
####### "PROVEAN_pred","M_CAP_pred","FATHMM_pred","SIFT4G_pred","SIFT_pred","PrimateAI_pred","DEOGEN2_pred")
####### scorevar<-c("MutPred_rankscore","REVEL_rankscore","VEST4_rankscore")
####### outfile<-paste0("/medpop/afib/data/NHLBI_WGS/freeze.8/Annotation/WGSA_parsed/s0/chr",num,"_dbNSFP_scatter_piece.gathered.sorted.select.tsv")



library(data.table)

dbnsfpfilterexport<-function(filename,normvar,scorevar=NULL,phredvar=NULL,outfile){

name1<-fread(cmd=paste0("zcat ",filename,"| head -n 2"),header=T,data.table=F,sep="\t")
name2<-names(name1)
colnumbers<-which(name2 %in% c(normvar,scorevar,phredvar))

com1<-paste0("zcat ",filename," | cut -f ",paste(colnumbers,collapse=",")," > ",outfile)
system(com1)
}



#######
####### create deleturious score
dbnsfpscore<-function(filename,normvar,scorevar=NULL,phredvar=NULL){

dat1<-dbnsfpfilter(filename=filename,normvar=normvar,scorevar=scorevar,phredvar=phredvar)

scorenewvars<-gsub("rankscore","pred",scorevar)
phrednewvars<-gsub("phred","pred",phredvar)

prevars<-c(normvar[grep("_pred",normvar)],scorenewvars,phrednewvars)

dat2<-dat1
for( ii in 1:length(prevars)){

prevar<-prevars[ii]
if(prevar %in% c("MetaLR_pred","MetaSVM_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","PROVEAN_pred","M_CAP_pred","FATHMM_pred","SIFT4G_pred","SIFT_pred","PrimateAI_pred","DEOGEN2_pred","MutPred_pred","REVEL_pred","VEST4_pred")){

dat2[,prevar]<-ifelse(is.na(dat2[,prevar]),NA,ifelse(dat2[,prevar]=="",NA,ifelse(dat2[,prevar]=="D",1,0)))

}

if(prevar %in% c("MutationAssessor_pred")){

dat2[,prevar]<-ifelse(is.na(dat2[,prevar]),NA,ifelse(dat2[,prevar]=="",NA,ifelse(dat2[,prevar]=="H",1,0)))

}

}

dat2$nmiss<-apply(dat2[,prevars],1,function(x){sum(!is.na(x))})
dat2$dscore<-apply(dat2[,prevars],1,function(x){sum(x,na.rm=T)})
dat2$dprop<-dat2$dscore/dat2$nmiss

return(dat2)
}
