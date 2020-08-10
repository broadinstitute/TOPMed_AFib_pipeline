####### no tissue weight * no prop deleterious score = weight = 1
####### prop cut off 0.1 - 1
opt.ptex.del.cutoff.SMMAT<-
function(num,gdsfile,groupfile,txannotfile,tissuename,aweight,phenfile,nullfile,geneid){

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
combphen2<-combphen[samples,]

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

##### plan 1
###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")

###### expression weight
nwei=1
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

#### del score weight
nwei=2
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-1;#txannot[,tissuename0]

###### new variable
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")

###### proportion cutoffs
prop0<-seq(0,1,by=0.1)
res0<-NULL

####### loop start
for (dw in 1:length(prop0)){
threprop<-prop0[dw]
txannot2<-subset(txannot,txannot[,tissuename0]>=threprop)
txannot2[,tissuedelweight]<-txannot2[,tissueweight[1]]*txannot2[,tissueweight[2]]
txannot2<-txannot2[,c("varid","ensg",tissuedelweight)]

######
###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot2,by.x=c("group_id","varid"),by.y=c("ensg","varid"))
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

###### annotation for TTN
optgene<-subset(annot3,group_id==geneid)

cat(threprop," starts\n")
######
###### grouping file
gr<-aggregateGRangesList(optgene)

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
res1$cuttoff<-threprop
cat(paste0(threprop," is done\n"))
print(res1)
res0<-rbind(res0,res1)
}

seqClose(seqData)
print(res0)
return(res0)
}


####### no tissue weight * no prop deleterious score = weight = 1
####### prop cut off 0.1 - 1
opt.del.cutoff.SMMAT<-
function(num,gdsfile,groupfile,txannotfile,tissuename,aweight,phenfile,nullfile,geneid){

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
combphen2<-combphen[samples,]

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

##### plan 1
###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")

###### expression weight
nwei=1
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

#### del score weight
nwei=2
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-txannot[,tissuename0]

###### new variable
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")

###### proportion cutoffs
prop0<-seq(0,1,by=0.1)
res0<-NULL

####### loop start
for (dw in 1:length(prop0)){
threprop<-prop0[dw]
#txannot2<-subset(txannot,txannot[,tissuename0]>=threprop)
#txanno2[,tissuedelweight]<-1
#txannot2<-txannot2[,c("varid","ensg",tissuedelweight)]

######
###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot,by.x=c("group_id","varid"),by.y=c("ensg","varid"),all=T)
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

###### annotation for TTN
optgene<-subset(annot3,group_id==geneid)

cat(threprop," starts\n")
######
###### grouping file
gr<-aggregateGRangesList(optgene)

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
res1$cuttoff<-threprop
cat(paste0(threprop," is done\n"))
print(res1)
res0<-rbind(res0,res1)
}

seqClose(seqData)
print(res0)
return(res0)
}



####### tissue weight * prop deleterious score
####### prop cut off 0.1 - 1
opt.del.noweight.SMMAT<-
function(num,gdsfile,groupfile,txannotfile,tissuename,aweight,phenfile,nullfile,geneid){

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
combphen2<-combphen[samples,]

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

##### plan 1
###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")

###### expression weight
nwei=1
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

#### del score weight
nwei=2
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-txannot[,tissuename0]

###### new variable
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")

###### proportion cutoffs
prop0<-seq(0,1,by=0.1)
res0<-NULL

####### loop start
for (dw in 1:length(prop0)){
threprop<-prop0[dw]
txannot2<-subset(txannot,txannot[,tissuename0]>=threprop)
txannot2[,tissuedelweight]<-txannot2[,tissueweight[1]]*txannot2[,tissueweight[2]]
txannot2<-txannot2[,c("varid","ensg",tissuedelweight)]

######
###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot2,by.x=c("group_id","varid"),by.y=c("ensg","varid"))
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

###### annotation for TTN
optgene<-subset(annot3,group_id==geneid)

cat(threprop," starts\n")
######
###### grouping file
gr<-aggregateGRangesList(optgene)

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
res1$cuttoff<-threprop
cat(paste0(threprop," is done\n"))
print(res1)
res0<-rbind(res0,res1)
}

seqClose(seqData)
print(res0)
return(res0)
}


####### tissue weight * prop deleterious score weight
####### alpha weight  1 - 10
opt.del.weight.SMMAT<-
function(num,gdsfile,groupfile,txannotfile,tissuename,aweight,phenfile,nullfile,geneid){

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
combphen2<-combphen[samples,]

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

##### plan 2

delweights<-c(1:10)
ares0<-NULL
###### start loops
for (dw in 1:length(delweights)){

print(dw)
aweight[2]<-delweights[dw]

###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")

###### tissue weight
nwei=1
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

#### del score weight
nwei=2
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
del_aweight<-aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

##### new variable
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")

##### tissue weight * delscore weight
txannot2<-txannot
txannot2[,tissuedelweight]<-txannot2[,tissueweight[1]]*txannot2[,tissueweight[2]]
txannot2<-txannot2[,c("varid","ensg",tissuedelweight)]


######
###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot2,by.x=c("group_id","varid"),by.y=c("ensg","varid"))
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

###### annotation for TTN
optgene<-subset(annot3,group_id==geneid)


cat(del_aweight," starts\n")

######
###### grouping file
gr<-aggregateGRangesList(optgene)

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
ares1<-assoc$results
ares1$alpha<-del_aweight
cat(paste0(del_aweight," is done\n"))
print(ares1)
ares0<-rbind(ares0,ares1)
}

seqClose(seqData)
print(ares0)
return(ares0)
}



####### tissue weight * prop deleterious score weight
####### tissue weight * prop deleterious score weight cut off 10% - 90%, 95%
opt.del.weight.prop.SMMAT<-
function(num,gdsfile,groupfile,txannotfile,tissuename,nweight,aweight,phenfile,nullfile,geneid){

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
combphen2<-combphen[samples,]

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

###### tissue type tissuename="Heart_Atrial_Appendage"
tissueweight<-paste0(tissuename,"_weight")

#### tissue weight
nwei=1
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

#### del prop weight
nwei=2
tissuename0<-tissuename[nwei]
tissueweight0<-tissueweight[nwei]
del_aweight<-aweight0<-aweight[nwei]
txannot[,tissuename0]<-as.numeric(txannot[,tissuename0])
txannot[,tissuename0]<-ifelse(txannot[,tissuename0]>1,NA,txannot[,tissuename0])
txannot[,tissueweight0]<-dbeta(txannot[,tissuename0],aweight0,1)

##### new variable
tissuedelweight<-paste0(c(tissuename,"weight"),collapse="_")

##### weight
txannot2<-txannot
txannot2[,tissuedelweight]<-txannot2[,tissueweight[1]]*txannot2[,tissueweight[2]]
txannot2<-txannot2[,c("varid","ensg",tissuedelweight)]

###### annotation file
annot<-get(load(groupfile))
annot$varid<-paste(annot$chr,annot$pos,annot$ref,annot$alt,sep=":")
annot2<-merge(annot,txannot2,by.x=c("group_id","varid"),by.y=c("ensg","varid"))
annot2<-annot2[,c(1,3:6,2,7:ncol(annot2))]
annot3<-subset(annot2,!is.na(annot2[,c(tissuedelweight)]))

###### annotation for TTN
optgene<-subset(annot3,group_id==geneid)

cutoffs<-quantile(optgene[,tissuedelweight],c(seq(0.1,0.9,0.1),0.95))

bres0<-NULL
for (dw in 1:length(cutoffs)){

cutoff<-cutoffs[[dw]]
cutperc<-names(cutoffs)[dw]

######
###### grouping file
optgene2<-subset(optgene,optgene[tissuedelweight]>=cutoff)
gr<-aggregateGRangesList(optgene2)

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
bres1<-assoc$results
bres1$perc<-cutperc
bres1$cutoffs<-cutoff
cat(paste0(cutperc," is done\n"))
print(bres1)
bres0<-rbind(bres0,bres1)
}

seqClose(seqData)
print(bres0)
return(bres0)
}
