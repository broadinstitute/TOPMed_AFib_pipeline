####### num=21
####### groupfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/var_grouping/hclof_missense/f8_hclof_missense_chr_",num,".RData")
####### gdsfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
####### varfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/varlist/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
####### phenfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/phenotype/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
####### nullfile<-"/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/nullmodel/AF_maleassocatedPCs_kinship.RData"
####### txannotfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/dbnsfp/TOPMed_Freeze8_LOF_missense_30delscore_chr",num,"_trinfo.RData")
####### tissuename<-c("Dprop")
####### cutoff<-0.9
####### geneid<-"ENSG00000134571"
####### outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/results/hclof_missense/Heart_Left_Ventricle/TOPMed_freeze8_AF_hclof_missense_noweight_",paste0(tissuename,collapse="_"),"_30delscore_noweighted_cutoff",cutoff,"_chr",num,".tsv")
####### results<-del.cutoff.SMMAT(num,gdsfile,groupfile,txannotfile,tissuename,cutoff,phenfile,nullfile)
####### write.table(results,outfile,col.names=T,row.names=F,quote=F,sep="\t")

del.cutoff.SMMAT<-
function(num,gdsfile,groupfile,txannotfile,tissuename,cutoff,phenfile,nullfile){

##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
names(phen1)[1]<-"sample.id"
samid0<-phen1$sample.id

######
###### QCed variants
vardata<-fread(varfile,header=F,sep="\t",data.table=F)
varid0<-vardata$V1

######
###### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)
samples <- seqGetData(gds, "sample.id")
missamples<-samples[!samples %in% samid0]
misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
colnames(misphen)<-names(phen1)
misphen$sample.id<-missamples
combphen<-rbind(phen1,misphen)
rownames(combphen)<-combphen$sample.id
combphen2<-combphen[as.character(samples),]

######
###### # construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

######
###### filter the gdsfile
seqSetFilter(seqData, sample.id=samid0, variant.id=varid0)

######
###### load null model
nullmod<-get(load(nullfile))

###### read transcript expression information
txannot<-get(load(txannotfile))

#### del score weight
nwei=1
tissuename0<-tissuename[nwei]
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissuedelweight]<-txannot[,tissuename0]

######
###### annotation file
threprop<-cutoff

annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot,by.x=c("group_id","varid"),by.y=c("ensg","varid"),all=T)
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
names(annot2)[1:5]<-c("group_id","chr","pos","ref","alt")
annot2[is.na(annot2$chr),"chr"]<-annot2[is.na(annot2$chr),"chr.y"]
annot2[is.na(annot2$pos),"pos"]<-annot2[is.na(annot2$pos),"pos.y"]
annot2[is.na(annot2$ref),"ref"]<-annot2[is.na(annot2$ref),"ref.y"]
annot2[is.na(annot2$alt),"alt"]<-annot2[is.na(annot2$alt),"alt.y"]
annot2$pos<-as.numeric(annot2$pos)
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

cat("starts\n")
######
###### grouping file
gr<-aggregateGRangesList(annot3)

######
###### create the iterator
iterator <- SeqVarListIterator(seqData, variantRanges=gr)

###### perfrom assocation test
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.001, test="SMMAT", weight.user=tissuedelweight, verbose=TRUE)

###### add one variable
assoc$results$cutoff<-threprop

######
cat(paste0(threprop," is done\n"))

######
seqClose(seqData)
return(assoc)
}
