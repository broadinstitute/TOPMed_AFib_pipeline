####
library(SeqArray)
library(SeqVarTools)
library(data.table)
library(GenomicRanges)
library(GENESIS)
library(GWASTools)
library(TopmedPipeline)


##### function start
hclofburden<-function(num,gdsfile,groupfile,phenfile,nullfile,outfile){

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
###### annotation file
annot<-get(load(groupfile))

######
###### grouping file
gr<-aggregateGRangesList(annot)

######
###### create the iterator
iterator <- SeqVarListIterator(seqData, variantRanges=gr)

######
###### load null model
nullmod<-get(load(nullfile))

###### perfrom assocation test
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.001, test="Burden", verbose=TRUE)

#####
save(assoc,file=outfile)
#####
seqClose(gds)
}


#######
####### summary of association results

summarydata<-function(files,chrs,thre_cMAC=0){
sumres<-NULL
for (num in 1:length(files)){
outfile<-files[num]
chr0<-chrs[num]
res1<-get(load(outfile))
sum0<-res1$results
sum0$chr<-chr0
varinfo0<-res1$variantInfo
cMAC0<-NULL
mpos0<-NULL
for (genenum in c(1:length(varinfo0))){
cMAC1<-sum(varinfo0[[genenum]]$MAC)
mpos1<-mean(varinfo0[[genenum]]$pos)
cMAC0<-c(cMAC0,cMAC1)
mpos0<-c(mpos0,mpos1)
}
sum0$cMAC<-cMAC0
sum0$mpos<-mpos0
sum1<-subset(sum0,cMAC>=thre_cMAC)
sumres<-rbind(sumres,sum1)
}
return(sumres)
}

#######
####### manhattan plot
#### manhattan plot for Jun 02 2018
man<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
    genomewideline = -log10(5e-08), highlight = NULL, highlight2=NULL, logp = TRUE,
    ...)
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x)))
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x)))
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x)))
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x)))
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]]))
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]]))
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]]))
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]]))
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        options(scipen = 999)
        d$pos = d$BP/1e+06
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index ==
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP +
                  lastbase
            }
            ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR ==
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.001)
    xmin = floor(max(d$pos) * -0.001)
#    def_args <- list(xaxt = "n",yaxt = "n", bty = "l",  xaxs = "i", yaxs = "i",
#        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
#            ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10]("p-value")),mgp=c(2,.5,0))
    def_args <- list(xaxt = "n",yaxt = "n", bty = "l",  xaxs = "i", yaxs = "i",
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
            ceiling(max(d$logp)+0.5)), xlab = "", ylab = "",mgp=c(2,.5,0))


    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1,cex.axis=0.8,tck=-0.03, line = 0, ...)
    }
    else {
        axis(1,cex.axis=0.8,tck=-0.03, line = 0, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
         for (i in unique(d$index)) {

                 with(d[d$index == unique(d$index)[i], ], points(pos,
                    logp, col = col[icol], pch = 20, ...))

     	#d1<-subset(d,logp<=20)
      #d2<-subset(d,logp>=45)
      #      with(d1[d1$index == unique(d1$index)[i], ], points(pos,
      #          logp, col = col[icol], pch = 20, ...))
      #      with(d2[d2$index == unique(d2$index)[i], ], points(pos,
      #          logp-25, col = col[icol], pch = 20, ...))
            icol = icol + 1

        }
    }
    if (suggestiveline)
        abline(h = suggestiveline, col = "blue",lty=3)
    if (genomewideline)
        abline(h = genomewideline, col = "black",lty=3)
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")

      d.highlight = d[which(d$SNP %in% highlight), ]
      with(d.highlight, points(pos, logp, col = "#377eb8", pch = 20, ...))

      #  d1<-subset(d,logp<=20)
      #  d2<-subset(d,logp>=45)
      #  d1.highlight = d1[which(d1$SNP %in% highlight), ]
      #  with(d1.highlight, points(pos, logp, col = "#377eb8", pch = 20, ...))
      #  d2.highlight = d2[which(d2$SNP %in% highlight), ]
      #  with(d2.highlight, points(pos, logp-25, col = "#377eb8", pch = 20,...))
    }
       if (!is.null(highlight2)) {
		if (any(!(highlight2 %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")

            d.highlight2 = d[which(d$SNP %in% highlight2), ]
            with(d.highlight2, points(pos, logp, col = "#e41a1c", pch = 20, ...))

      #  d1<-subset(d,logp<20)
		  # d1.highlight2 = d1[which(d1$SNP %in% highlight2), ]
      #  with(d1.highlight2, points(pos, logp, col="#e41a1c", pch = 20, ...))

    }

}


#####
##### allele count calculation

countcal<-function(num,gdsfile,phenfile,outfile){

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

######
###### filter the gdsfile
seqSetFilter(gds, sample.id=samid0, variant.id=varid0)

######
###### total counts
count0<-SeqVarTools::alleleCount(gds, n=0, use.names=FALSE)
count1<-SeqVarTools::alleleCount(gds, n=1, use.names=FALSE)

######
###### case samples
case0<-subset(phen1,AF_status=="case")
caseid0<-case0$sample.id

###### case counts
seqSetFilter(gds, sample.id=caseid0, variant.id=varid0)
case.count0<-SeqVarTools::alleleCount(gds, n=0, use.names=FALSE)
case.count1<-SeqVarTools::alleleCount(gds, n=1, use.names=FALSE)

######
###### control samples
cont0<-subset(phen1,AF_status=="control")
contid0<-cont0$sample.id

###### case counts
seqSetFilter(gds, sample.id=contid0, variant.id=varid0)
cont0.count0<-SeqVarTools::alleleCount(gds, n=0, use.names=FALSE)
cont0.count1<-SeqVarTools::alleleCount(gds, n=1, use.names=FALSE)

###### count data
countdata<-data.frame(variant.id=varid0,total.ref=count0,total.alt=count1,case.ref=case.count0,case.alt=case.count1,cont.ref=cont0.count0,cont.alt=cont0.count1)

######
###### write table
write.table(countdata,outfile,col.names=T,row.names=F,quote=F,sep="\t")

#####
#####
sessionInfo()
quit("no")

}


####
#### single variant association tests
singleassoc<-function(num=num,gdsfile=gdsfile,phenfile=phenfile,nullfile=nullfile,cutoff=0.001,outfile=outfile){

##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
names(phen1)[1]<-"sample.id"
samid0<-phen1$sample.id

######
###### QCed variants
vardata<-fread(varfile,header=T,sep="\t",data.table=F,select=c(1:3))
vardata$freq<-vardata$total.alt/(vardata$total.ref+vardata$total.alt)
vardata2<-subset(vardata,freq >=cutoff & freq<=(1-cutoff))
varid0<-vardata2$variant.id

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
###### create the iterator
iterator <- SeqVarBlockIterator(seqData, verbose=TRUE)

######
###### load null model
nullmod<-get(load(nullfile))

######
###### perfrom test
assoc <- assocTestSingle(iterator, nullmod, test="Score",verbose=TRUE)

######
###### save files
write.table(assoc,outfile,col.names=T,row.names=F,quote=F,sep="\t")

#####
#####
sessionInfo()
quit("no")

}


####
#### rare variant testing for variance

omnibus_test<-function(vartest,num,gdsfile,groupfile,phenfile,nullfile,outfile){

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
###### annotation file
annot<-get(load(groupfile))

######
###### grouping file
gr<-aggregateGRangesList(annot)

######
###### create the iterator
iterator <- SeqVarListIterator(seqData, variantRanges=gr)

######
###### load null model
nullmod<-get(load(nullfile))

###### perfrom assocation test
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.001, weight.beta = c(1,25), test=vartest, verbose=TRUE)

#####
save(assoc,file=outfile)
#####
seqClose(gds)
}

#####
##### carrier information

carrier_info<-function(phenfile,num,gdsfile,geneid,resultfile){

##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
names(phen1)[1]<-"sample.id"
samid0<-phen1$sample.id

######
###### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)

######
###### read result file
res<-get(load(resultfile))
varinfo<-res$variantInfo[[geneid]]
varid0<-varinfo$variant.id

#####
##### filter
seqSetFilter(gds, sample.id=samid0, variant.id=varid0)

####
#### dosage
dosage0<-altDosage(gds)

####
#### carriers
carrierlist<-apply(dosage0,2,function(x){names(which(x>=1))})
varids<-names(carrierlist)
comb<-NULL
for ( ii in 1:length(carrierlist)){
ids<-carrierlist[[ii]]
sub1<-subset(phen1,sample.id %in% ids)
sub1$varid<-varids[ii]
comb<-rbind(comb,sub1)
}
seqClose(gds)
return(comb)
}
