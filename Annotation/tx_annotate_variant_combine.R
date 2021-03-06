### example
R.utils::sourceDirectory("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/TOPMed_AFib_pipeline/R/", modifiedOnly=TRUE);


####### tx annotataion
####### hc_lof and missense combine
loffile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx/TOPMed_Freeze8_LOF_chrall_trinfo.tsv.gz"
missesefile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx/TOPMed_Freeze8_missense_chrall_trinfo.tsv.gz"
outfiles<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx/TOPMed_Freeze8_LOF_missense_chr",c(1:22,"X"),"_trinfo.RData")

####### step 1
lof0<-fread(cmd=paste0("zcat ",loffile),header=T,data.table=F,sep="\t")
missese0<-fread(cmd=paste0("zcat ",missesefile),header=T,data.table=F,sep="\t")


for (num in c(2:22,"X")){
chrnum<-paste0("chr",num)
outfile<-outfiles[as.numeric(num)]
lof1<-subset(lof0,CHROM==chrnum)
missese1<-subset(missese0,CHROM==chrnum)

comb<-rbind(lof1,missese1)
comb1<-comb[order(comb$POSITION),]
save(comb1,file=outfile)
}


#################
################# deleterious score and tx Annotation
num=1
for (num in c(1:22,"X")){

scorefile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/s0/chr",num,"_dbNSFP_scatter_piece.gathered.sorted.pred.tsv.gz")
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx/TOPMed_Freeze8_LOF_missense_chr",num,"_trinfo.RData")
txdelannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")

###### score file
library(tidyr)
score1<-fread(cmd=paste0("gunzip -c ",scorefile),header=T,data.table=F,sep="\t")
score1$vargid<-paste(score1$chr,score1$pos,score1$ref,score1$alt,score1$Ensembl_geneid,sep=":")
score2<-subset(score1,nmiss>=5)
score3<-score2[order(score2$vargid,-score2$dprop),]
score4<-score3[!duplicated(score3$vargid),]

#######
tx0<-get(load(txannotfile))
tx0$vargid<-paste0(tx0$varid,":",tx0$ensg)
tx_score<-merge(tx0,score4,by="vargid",all=T)

#######
tx_score<-separate(data=tx_score,col="vargid",c("chr","pos","ref","alt","ensg"),remove=F)
tx_score[which(tx_score$lof=="HC"),"dprop"]<-1

write.table(tx_score,txdelannotfile,col.names=T,row.names=F,quote=F,sep="\t")
}

##### TTN
##### input files
num=2
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000155657"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_noweighted_cutoff_0_1_gene_",geneid,".tsv")

########
########
result0<-opt.del.noweight.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result0,outfile,col.names=T,row.names=F,quote=F,sep="\t")

##### LMNA
##### input files
num=1
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000160789"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_noweighted_cutoff_0_1_gene_",geneid,".tsv")

########
########
result0<-opt.del.noweight.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result0,outfile,col.names=T,row.names=F,quote=F,sep="\t")

##### MYBPC3
##### input files
num=11
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000134571"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_noweighted_cutoff_0_1_gene_",geneid,".tsv")

########
########
result0<-opt.del.noweight.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result0,outfile,col.names=T,row.names=F,quote=F,sep="\t")


######### PLAN 1 result figure
geneids<-c("ENSG00000155657","ENSG00000160789","ENSG00000134571")
outfiles<-paste("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_noweighted_cutoff_0_1_gene_",geneids,".tsv")

outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_noweighted_cutoff_0_1_gene_",geneid,".tsv")




############
############ PLAN2 wighted on del score

##### TTN
##### input files
num=2
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000155657"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_weighted_alpha_1_10_gene_",geneid,".tsv")

result1<-opt.del.weight.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result1,outfile,col.names=T,row.names=F,quote=F,sep="\t")



##### LMNA
##### input files
num=1
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000160789"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_weighted_alpha_1_10_gene_",geneid,".tsv")


result1<-opt.del.weight.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result1,outfile,col.names=T,row.names=F,quote=F,sep="\t")


##### MYBPC3
##### input files
num=11
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000134571"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_weighted_alpha_1_10_gene_",geneid,".tsv")


result1<-opt.del.weight.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result1,outfile,col.names=T,row.names=F,quote=F,sep="\t")

########
##### input files TTN
num=2
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7,7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000155657"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_weighted_",aweight[2],"cutoff_10p_90p_95p_gene_",geneid,".tsv")

result2<-opt.del.weight.prop.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result2,outfile,col.names=T,row.names=F,quote=F,sep="\t")


#######
####### LMNA
num=1
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7,7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000160789"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_weighted_",aweight[2],"cutoff_10p_90p_95p_gene_",geneid,".tsv")

result2<-opt.del.weight.prop.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result2,outfile,col.names=T,row.names=F,quote=F,sep="\t")

######
###### MYBPC3
num=11
groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/tx_del/TOPMed_Freeze8_LOF_missense_delscore_chr",num,"_trinfo.RData")
aweight<-c(7,7)
tissuename<-c("Heart_Left_Ventricle","dprop")
geneid<-"ENSG00000134571"
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_",paste0(tissuename,collapse="_"),"_delscore_weighted_",aweight[2],"cutoff_10p_90p_95p_gene_",geneid,".tsv")

result2<-opt.del.weight.prop.SMMAT(num=num,gdsfile=gdsfile,groupfile=groupfile,txannotfile=txannotfile,tissuename=tissuename,aweight=aweight,phenfile=phenfile,nullfile=nullfile,geneid=geneid)
write.table(result2,outfile,col.names=T,row.names=F,quote=F,sep="\t")






#res0<-NULL
delweights<-c(1:10)
#for (dw in 1:length(delweights)){
res0<-foreach(dw=1:length(delweights), .combine=rbind) %dopar% {

print(dw)
aweight[2]<-delweights[dw]

###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")
for (nwei in 1:length(tissuename)){
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)
}
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")
txannot[,tissuedelweight]<-txannot[,tissueweight[1]]*txannot[,tissueweight[2]]
txannot2<-txannot[,c("varid","ensg",tissuedelweight)]
txannot2<-subset(txannot2,txannot2[,tissuedelweight]>=0.2)

######
###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot2,by.x=c("group_id","varid"),by.y=c("ensg","varid"))
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

######
###### grouping file
gr<-aggregateGRangesList(annot3)

######
###### create the iterator
iterator <- SeqVarListIterator(seqData, variantRanges=gr)


###### perfrom assocation test
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.001, test="SMMAT", weight.user=tissuedelweight, verbose=TRUE)

######
###### # construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

######
###### filter the gdsfile
seqSetFilter(seqData, sample.id=samid0, variant.id=varid0)

##### results
res1<-assoc$results

res1
}
res0



#####
save(assoc,file=outfile)
#####
seqClose(gds)
}
>
####
git clone https://github.com/UW-GAC/analysis_pipeline.git
