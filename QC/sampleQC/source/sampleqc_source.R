##### setup library
library(SeqVarTools)
library(parallel)
library(data.table)

.nSamp <- function(gdsobj) {
  #sum(seqGetFilter(gdsobj)$sample.sel)
  seqSummary(gdsobj, "genotype", check="none", verbose=FALSE)$seldim[2L]
}

##### heterozygosity counts
parahetcount <- function(gdsobj,parallel=seqGetParallel()){

het <- integer(.nSamp(gdsobj))
nonmiss <- integer(.nSamp(gdsobj))

seqParallel(parallel, gdsobj, FUN = function(f) {
seqApply(f, "genotype", function(x) {
                       nm <- !is.na(x[1,]) & !is.na(x[2,])
                       het <<- het + (x[1,] != x[2,] & nm)
                       nonmiss <<- nonmiss + nm
                     },
                   margin="by.variant", as.is="none")})
het<-data.frame(het, nonmiss)
return(het)
}

##### non-reference homozygosity counts
parahomcount<-function(gdsobj, allele=c("any", "ref", "alt"),
           margin=c("by.variant", "by.sample"), parallel=seqGetParallel(),use.names=FALSE) {
    hom.func <- switch(match.arg(allele),
                       any=function(a,b) {a == b},
                       ref=function(a,b) {a == b & a == 0},
                       alt=function(a,b) {a == b & a > 0})
    margin <- match.arg(margin)
    if (margin == "by.variant") {
      hom <- seqParallel(parallel, gdsobj, FUN = function(f) {
      seqApply(f, "genotype",function(x) {
                        sum(hom.func(x[1,], x[2,]), na.rm=TRUE) /
                          sum(!is.na(x[1,]) & !is.na(x[2,]))
                        },
                      margin=margin,as.is="double")})
      if (use.names)
        names(hom) <- seqGetData(gdsobj, "variant.id")
      hom
    } else {
      hom <- integer(.nSamp(gdsobj))
      nonmiss <- integer(.nSamp(gdsobj))
      ## use "<<-" operator to find "hom" in the parent environment
      seqParallel(parallel, gdsobj, FUN = function(f) {
      seqApply(f, "genotype",function(x) {
                 nm <- !is.na(x[1,]) & !is.na(x[2,])
                 hom <<- hom + (hom.func(x[1,], x[2,]) & nm)
                 nonmiss <<- nonmiss + nm
               },
             margin="by.variant", as.is="none")})
      if (use.names)
        names(hom) <- seqGetData(gdsobj, "sample.id")
        homout<-data.frame(hom, nonmiss)
    }
  }


  .isTransition <- function(ref, alt) {
    (ref %in% c("C","T") & alt %in% c("C","T")) |
    (ref %in% c("A","G") & alt %in% c("A","G"))
  }

  .isTransversion <- function(ref, alt) {
    (ref %in% c("C","T") & alt %in% c("A","G")) |
    (ref %in% c("A","G") & alt %in% c("C","T"))
  }




  paratitvcount<-function(gdsobj, by.sample=FALSE, parallel=seqGetParallel(),use.names=FALSE) {
    ref <- refChar(gdsobj)
    alt <- altChar(gdsobj)
    ti <- .isTransition(ref, alt)
    tv <- .isTransversion(ref, alt)
    if (by.sample) {
      isVar <- function(x) {colSums(x, na.rm=TRUE) > 0}
      tisum <- integer(.nSamp(gdsobj))
      tvsum <- integer(.nSamp(gdsobj))
      seqParallel(parallel, gdsobj, FUN = function(f) {
      seqApply(f, "genotype",
               function(index, x) {
                 tisum <<- tisum + (ti[index] & isVar(x));
                 tvsum <<- tvsum + (tv[index] & isVar(x))
               },
             margin="by.variant",as.is="none",var.index="relative")})
      titv <- data.frame(tisum,tvsum)
      if (use.names)
        row.names(titv) <- seqGetData(gdsobj, "sample.id")
      titv
    } else {
      sum(ti) / sum(tv)
    }
  }

####
####
sample_qc_chunk<-function(gdsfile,varfile,phenfile,chunknum=1,nvariants=50000){

cpunum<-1
#####
##### samples
phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
samid0<-phen1$SAMPLE_ID

######
###### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)

######
###### non non-nomomorphic variants
######
varlist0<-fread(varfile,header=F,sep="\t",data.table=F)
varlist0<-vardata$V1
nvar0<-length(varlist0)
varlist1<-c(seq(1,nvar0,by=nvariants),nvar0+1)
varlist<-varlist0[c((varlist1[chunknum]):(varlist1[chunknum+1]-1))]
varid0<-varlist

######
###### filter gds
seqSetFilter(gds, sample.id=samid0, variant.id=varid0)

######
###### setup the cluster
seqParallelSetup()

##### call rate
missingrate0<-seqMissing(gds, per.variant=FALSE, .progress=TRUE,parallel=cpunum)
names(missingrate0)<-samid0

##### SNP  Indel ratio
snp0<-isSNV(gds)
names(snp0)<-varid0
snpvar<-names(which(snp0==TRUE))
indelvar<-names(which(snp0==FALSE))

##### SNP
seqSetFilter(gds, variant.id=snpvar,sample.id=samid0)
snphet<-parahetcount(gds,parallel=cpunum)
snphom<-parahomcount(gds,allele="alt",margin="by.sample",parallel=cpunum)
nalt<-snphet$het+snphom$hom*2
snpout<-data.frame(snpalt=nalt,nsnp=snphet$nonmiss)
seqResetFilter(gds)

###### INDEL
seqSetFilter(gds, variant.id=indelvar,sample.id=samid0)
indelhet<-parahetcount(gds,parallel=cpunum)
indelhom<-parahomcount(gds,allele="alt",margin="by.sample",parallel=cpunum)
nalt<-indelhet$het+indelhom$hom*2
indelout<-data.frame(indelalt=nalt,nindel=indelhet$nonmiss)
seqResetFilter(gds)

##### HET nonrefhom count
seqSetFilter(gds, sample.id=samid0, variant.id=varid0)
allhet<-snphet$het+indelhet$het
allhom<-snphom$hom+indelhom$hom
nvar<-snphet$nonmiss+indelhet$nonmiss
hethomout<-data.frame(allhet,allhom,nvar)

##### TI TV ratio
titvcount<-paratitvcount(gds, by.sample=TRUE,parallel=cpunum,use.names=TRUE)

###### singleton
ac1count<-countSingletons(gds,use.names=TRUE)

##### total variants included in a chromosome
totalvar<-length(varid0)

###### create data set
out1<-data.frame(NWDID=samid0,missing=missingrate0,snpout,indelout,hethomout,titvcount,ac1count,totalvar)

#### session Info
seqClose(gds)
seqParallelSetup(FALSE)
sessionInfo()

return(out1)
}



####
#### chunknum
totalchunk<-function(varfile,nvariants=50000){
vardata<-fread(varfile,header=F,sep="\t",data.table=F)
varlist0<-vardata$V1
nvar0<-length(varlist0)
varlist1<-c(seq(1,nvar0,by=nvariants),nvar0+1)
totalchunk<-length(varlist1)-1
return(totalchunk)
}
