######## logitf
library(logistf)
library(data.table)
ageatonset<-function(phenfile,num,geneid,gdsfile,resultfile,unrelcol,noHF_col){

cars<-carrier_info(phenfile=phenfile,num=num,gdsfile=gdsfile,geneid=geneid,resultfile=resultfile)

######
###### HF free unrelated
pheno<-fread(phenfile,header=T,data.table=F,sep="\t")
pheno$carriers<-0
pheno[pheno$SAMPLE_ID %in% cars$sample.id,"carriers"]<-1
pheno$AF<-ifelse(pheno$AF_status=="case",1,0)
table(pheno$carriers)
pheno$carriers<-as.factor(pheno$carriers)

######
###### unrel
unrel<-subset(pheno,pheno[,unrelcol]==1 )
table(unrel$carriers)
######
######
cont<-subset(unrel,AF==0)
table(cont$carriers)

######
######
case0<-subset(unrel,AF==1)
table(case0$carriers)

case0$age_at_onset_new<-ifelse(is.na(case0$Age_at_AF_onset),100,case0$Age_at_AF_onset)

agecutoffs<-c(300,65,55,45,35)

####### unrelated
preres0<-NULL
for (ii in 1:length(agecutoffs)){

agecutoff<-agecutoffs[ii]
case1<-subset(case0,age_at_onset_new<agecutoff)
tab1<-table(case1$carriers)
ncarr0<-tab1[["1"]]
totaln<-sum(tab1)
prev1<-prop.table(tab1)[["1"]]
preres1<-data.frame(agecutoff,totaln,ncarr0,prev1,stringsAsFactors=F)
preres0<-rbind(preres0,preres1)
}

tab1<-table(cont$carriers)
ncarr0<-tab1[["1"]]
totaln<-sum(tab1)
prev1<-prop.table(tab1)[["1"]]
preres1<-data.frame(agecutoff="control",totaln,ncarr0,prev1,stringsAsFactors=F)
preres0<-rbind(preres0,preres1)

preres0$type<-"unrel"
unrelres<-preres0

#####
##### unrelated + no HF before AF

case01<-subset(case0,case0[,noHF_col]==0)
table(case01$carriers)

rmlist<-subset(case01,case01[,"LVEF_by_echo"]<50) ### nobody


####### unrelated + no HF before AF
preres0<-NULL
for (ii in 1:length(agecutoffs)){

agecutoff<-agecutoffs[ii]
case1<-subset(case01,age_at_onset_new<agecutoff)
tab1<-table(case1$carriers)
ncarr0<-tab1[["1"]]
totaln<-sum(tab1)
prev1<-prop.table(tab1)[["1"]]
preres1<-data.frame(agecutoff,totaln,ncarr0,prev1,stringsAsFactors=F)
preres0<-rbind(preres0,preres1)
}

tab1<-table(cont$carriers)
ncarr0<-tab1[["1"]]
totaln<-sum(tab1)
prev1<-prop.table(tab1)[["1"]]
preres1<-data.frame(agecutoff="control",totaln,ncarr0,prev1,stringsAsFactors=F)
preres0<-rbind(preres0,preres1)

preres0$type<-"unrel_noHF_before_AF"
unrel_noHF<-preres0
result0<-rbind(unrelres,unrel_noHF)

result0$geneid<-geneid

return(result0)
}
