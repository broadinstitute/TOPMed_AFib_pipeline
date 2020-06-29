protein_coding_variants<-function(num=num,gdsfile=gdsfile,phenfile=phenfile,outfile=outfile){

##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
names(phen1)[1]<-"sample.id"
samid0<-phen1$sample.id

###### variants
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
combphen2<-combphen[samples,]

######
###### # construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))

######
###### filter the gdsfile
seqSetFilter(seqData, sample.id=samid0, variant.id=varid0)

######
###### protien coding variants
ptregion<-get(load(file="/medpop/afib/schoi/projects/TOPMed/Freeze6/Data/ECG/annotation/Ensembl/ENSEMBL_exon_info_mart_export_Apr302019_protein_coding_region.RData"))
ptregion2<-subset(ptregion,ptregion$exonstart != ptregion$exonend)
ptregion2$exonstart<-ptregion2$exonstart-2
ptregion2$exonend<-ptregion2$exonend+2
ptregion3<-subset(ptregion2,chromosome_name==num)

if(num!="X"){
seqSetFilterChrom(seqData,include=as.numeric(ptregion3$chromosome_name),is.num=TRUE,from.bp=ptregion3$exonstart,to.bp=ptregion3$exonend,intersect=TRUE)
}else{
seqSetFilterChrom(seqData,include=as.character(ptregion3$chromosome_name),is.num=FALSE,from.bp=ptregion3$exonstart,to.bp=ptregion3$exonend,intersect=TRUE)
}

varids <- seqGetData(seqData, "variant.id")

write.table(varids,outfile,col.names=F,row.names=F,quote=F,sep="\t")
seqClose(gds)
}
