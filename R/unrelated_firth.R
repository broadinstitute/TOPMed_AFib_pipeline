######## logitf
library(logistf)
library(data.table)
unrelated_firth<-function(phenfile,num,geneid,gdsfile,resultfile){

cars<-carrier_info(phenfile=phenfile,num=num,gdsfile=gdsfile,geneid=geneid,resultfile=resultfile)

######
###### HF free unrelated
pheno<-fread(phenfile,header=T,data.table=F,sep="\t")
pheno$carriers<-0
pheno[pheno$SAMPLE_ID %in% cars$sample.id,"carriers"]<-1
pheno$AF<-ifelse(pheno$AF_status=="case",1,0)

######
###### unrel
unrel<-subset(pheno,unrel_3degree==1 )
unrel<-subset(unrel,HF_before_AF==0)

covs<-c("Sex","PC1","PC2","PC4","PC5","PC6","PC8","PC10","PC11")
form1<-formula(paste0("AF ~ carriers +",paste0(covs,collapse="+")))
mod1<-logistf(form1,data=unrel)

return(mod1)
}
