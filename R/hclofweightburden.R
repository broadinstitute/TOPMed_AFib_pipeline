####
#### LOF carrier status
hclofweightburden<-function(num,gdsfile,groupfile,txannotfile,tissuename,aweight,phenfile,nullfile,outfile){

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

###### read transcript expression information
txannot<-fread(txannotfile,header=T,data.table=F,sep="\t")

###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")
txannot[,tissuename]<-ifelse(txannot[,tissuename]>1,NA,txannot[,tissuename])
txannot[,tissueweight]<-dbeta(txannot[,tissuename],aweight,1)
txannot2<-txannot[,c("varid","ensg",tissueweight)]

######
###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot2,by.x=c("group_id","varid"),by.y=c("ensg","varid"))
annot2<-annot2[,c(1,3:6,2,7)]
annot3<-subset(annot2,!is.na(annot2[,c(tissueweight)]))
#annot3<-as_tibble(annot3)

######
###### grouping file
gr<-aggregateGRangesList(annot3)

######
###### create the iterator
iterator <- SeqVarListIterator(seqData, variantRanges=gr)

######
###### load null model
nullmod<-get(load(nullfile))

###### perfrom assocation test
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.001, test="Burden", weight.user=tissueweight, verbose=TRUE)

#####
save(assoc,file=outfile)
#####
seqClose(gds)
}
