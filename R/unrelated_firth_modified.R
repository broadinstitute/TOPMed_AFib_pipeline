######## logitf
library(logistf)
library(data.table)

unrelated_firth_standard<-function(phenfile,num,geneid,gdsfile,resultfile,disease="AF_status",cov=NULL,PC=NULL,unrelcol=NULL,standard=TRUE){
disease0<-disease

#### identify the carriers
cars<-carrier_info(phenfile=phenfile,num=num,gdsfile=gdsfile,geneid=geneid,resultfile=resultfile)

######
###### HF free unrelated
pheno<-fread(phenfile,header=T,data.table=F,sep="\t")
names(pheno)[1]<-"sample.id"
pheno$carriers<-0
pheno[pheno$sample.id %in% cars$sample.id,"carriers"]<-1
print(nrow(pheno))

#####
if(class(pheno[,disease0])=="character"){
pheno$disease<-ifelse(pheno[,disease0]=="case",1,0)
}else{
pheno$disease<-pheno[,disease0]
}

######
###### unrelaed samples
if(!is.null(unrelcol){
pheno<-subset(pheno,pheno[,unrelcol]==1)
}
print(nrow(pheno))

######
###### standardize the PCs
if(length(PC)>0){
if(standard){
######
###### pcs will be standardized
pheno[,PC]<-apply(pheno[,PC],2,function(x){(x-mean(x))/sd(x)})
}else{}
}else{}

covs<-c(cov,PC)
form1<-formula(paste0("disease ~ carriers +",paste0(covs,collapse="+")))

###### unrelated samples
mod1<-logistf(form1,data=pheno)
mod1$data<-NULL

##### save result
res<-list(unrelated=mod1)
return(res)
}
