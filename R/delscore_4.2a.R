
delscore4.2a<-function(filename,delvar,mintool=7){
library(data.table)
library(tidyr)

dat0<-fread(cmd=paste0("zcat ",filename),header=T,data.table=F,sep="\t")
vars<-c("MetaRNN_pred","SIFT_pred","SIFT4G_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred","MutationAssessor_pred","FATHMM_pred","PROVEAN_pred","VEST4_rankscore","MetaSVM_pred","MetaLR_pred","M-CAP_pred","REVEL_rankscore","MutPred_rankscore","MVP_rankscore","MPC_rankscore","PrimateAI_pred","DEOGEN2_pred","BayesDel_addAF_pred","BayesDel_noAF_pred","ClinPred_pred","LIST-S2_pred","Aloft_pred","Aloft_Confidence","CADD_phred","DANN_rankscore","fathmm-MKL_coding_pred","fathmm-XF_coding_pred","Eigen-phred_coding","Eigen-PC-phred_coding")
varlist<-names(dat0)
vars<-varlist[varlist %in% vars]

##########
##########
dselect<-c("MetaRNN_pred","SIFT_pred","SIFT4G_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred","FATHMM_pred","PROVEAN_pred","MetaSVM_pred","MetaLR_pred","M-CAP_pred","PrimateAI_pred","DEOGEN2_pred","BayesDel_addAF_pred","BayesDel_noAF_pred","ClinPred_pred","LIST-S2_pred","fathmm-MKL_coding_pred","fathmm-XF_coding_pred")
dselect<-vars[vars %in% dselect]
for (ii in 1:length(dselect)){
variable0<-dselect[ii]
table(dat0[,variable0])
newvar<-paste0(variable0,"_D")
dat0[,newvar]<-NA
dat0[grep("T|P|B|N|U",dat0[,variable0]),newvar]<-0
dat0[grep("D",dat0[,variable0],fixed=T),newvar]<-1
##print(table(dat0[,newvar],useNA="always"))
}
##########
##########
hselect<-c("MutationAssessor_pred")
hselect<-vars[vars %in% hselect]
for (ii in 1:length(hselect)){
variable0<-hselect[ii]
##table(dat0[,variable0])
newvar<-paste0(variable0,"_D")
dat0[,newvar]<-NA
dat0[grep("M|L|N",dat0[,variable0]),newvar]<-0
dat0[grep("H",dat0[,variable0],fixed=T),newvar]<-1
##print(table(dat0[,newvar],useNA="always"))
}

#########
######### ranks
rank<-c("VEST4_rankscore","REVEL_rankscore","MutPred_rankscore","MVP_rankscore","MPC_rankscore","DANN_rankscore")
rank<-vars[vars %in% rank]

for (ii in 1:length(rank)){
variable0<-rank[ii]
newvar<-paste0(variable0,"_D")
dat0[,newvar]<-NA
dat0[,variable0]<-as.numeric(dat0[,variable0])
dat0[which(dat0[,variable0]<0.9),newvar]<-0
dat0[which(dat0[,variable0]>=0.9),newvar]<-1
#print(table(dat0[,newvar],useNA="always"))
}

########
######## phred scale
phred<-c("CADD_phred","Eigen-phred_coding","Eigen-PC-phred_coding")
phred<-vars[vars %in% phred]

for (ii in 1:length(phred)){
variable0<-phred[ii]
#summary(dat0[,variable0])
newvar<-paste0(variable0,"_D")
dat0[,newvar]<-NA
dat0[,variable0]<-as.numeric(dat0[,variable0])
#summary(dat0[,variable0])
dat0[which(dat0[,variable0]<10),newvar]<-0
dat0[which(dat0[,variable0]>=10),newvar]<-1
#print(table(dat0[,newvar],useNA="always"))
}

########
######## aloft is difficult
aloft<-c("Aloft_pred","Aloft_Confidence")
variable0<-aloft[1]
#table(dat0[,variable0])

newvar<-paste0(variable0,"_D")
dat0$taloft<-gsub("Tolerant","T",gsub("Recessive","D",gsub("Dominant","D",dat0[,aloft[1]])))
dat0$taloftconf<-gsub("Low Confidence","L",gsub("High Confidence","D",dat0[,aloft[2]]))
dat0[,newvar]<-NA
dat0[grep("T|D",dat0[,variable0]),newvar]<-0

test<-dat0[grep("D",dat0$taloft),]
for (ii in 1:nrow(test)){
trow<-test[ii,]
dlocdel<-grep("D",unlist(strsplit(trow$taloft,",")))
dlocconf<-grep("D",unlist(strsplit(trow$taloftconf,",")))
test[ii,newvar]<-ifelse(sum(dlocdel %in% dlocconf)>=1,1,test[ii,newvar])
}

dat0[rownames(test),newvar]<-test[,newvar]
#print(table(dat0[rownames(test),"taloftconf"],dat0[rownames(test),newvar]))
#print(table(dat0[,newvar],useNA="always"))

###########
########### create score
dvars<-paste0(vars[!vars %in% "Aloft_Confidence"],"_D")
dat1<-dat0[,dvars]
dat1$Dscore<-apply(dat1[,dvars],1,sum,na.rm=T)
dat1$Dnmiss<-apply(dat1[,dvars],1,function(x){sum(!is.na(x))})
dat1$Dprop<-apply(dat1[,dvars],1,function(x){mean(x,na.rm=T)})

###########
########### add the score to data
dat0$Dscore<-dat1$Dscore
dat0$Dtools<-dat1$Dnmiss
dat0$Dprop<-dat1$Dprop

##########
dat2<-subset(dat0,Dtools>=mintool)
dat2$gvarID<-paste0(dat2[,1],":",dat2[,4])
dat3<-dat2[order(dat2$gvarID,-dat2$Dprop),]
dat4<-dat3[!duplicated(dat3$gvarID),]
dat5<-dat4[,c("gvarID","Dscore","Dtools","Dprop")]
########## var deleterious info
names(dat5)[4]<-delvar

dat5<-separate(data=dat5,col="gvarID",into=c("chr","pos","ref","alt","geneID"))
#write.table(dat5,outfile,col.names=T,row.names=F,quote=F,sep="\t")
return(dat5)
}
