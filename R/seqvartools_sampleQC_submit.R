
.nSamp <- function(gdsobj) {
  #sum(seqGetFilter(gdsobj)$sample.sel)
  seqSummary(gdsobj, "genotype", check="none", verbose=FALSE)$seldim[2L]
}

##### heterozygosity counts
parahetcount <- function(gdsobj){

het <- integer(.nSamp(gdsobj))
nonmiss <- integer(.nSamp(gdsobj))

seqParallel(, gdsobj, FUN = function(f) {
seqApply(f, "genotype", as.is="none",
                     function(x) {
                       nm <- !is.na(x[1,]) & !is.na(x[2,])
                       het <<- het + (x[1,] != x[2,] & nm)
                       nonmiss <<- nonmiss + nm
                     })},split="by.variant")
het<-data.frame(het, nonmiss)
return(het)
}

##### non-reference homozygosity counts
parahomcount<-function(gdsobj, allele=c("any", "ref", "alt"),
           margin=c("by.variant", "by.sample"), use.names=FALSE) {
    hom.func <- switch(match.arg(allele),
                       any=function(a,b) {a == b},
                       ref=function(a,b) {a == b & a == 0},
                       alt=function(a,b) {a == b & a > 0})
    margin <- match.arg(margin)
    if (margin == "by.variant") {
      hom <- seqParallel(, gdsobj, FUN = function(f) {
      seqApply(f, "genotype",as.is="double",
                      function(x) {
                        sum(hom.func(x[1,], x[2,]), na.rm=TRUE) /
                          sum(!is.na(x[1,]) & !is.na(x[2,]))
                      })},split=margin)
      if (use.names)
        names(hom) <- seqGetData(gdsobj, "variant.id")
      hom
    } else {
      hom <- integer(.nSamp(gdsobj))
      nonmiss <- integer(.nSamp(gdsobj))
      ## use "<<-" operator to find "hom" in the parent environment
      seqParallel(, gdsobj, FUN = function(f) {
      seqApply(f, "genotype",as.is="none",
               function(x) {
                 nm <- !is.na(x[1,]) & !is.na(x[2,])
                 hom <<- hom + (hom.func(x[1,], x[2,]) & nm)
                 nonmiss <<- nonmiss + nm
               })},split="by.variant")
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




  paratitvcount<-function(gdsobj, by.sample=FALSE, use.names=FALSE) {
    ref <- refChar(gdsobj)
    alt <- altChar(gdsobj)
    ti <- .isTransition(ref, alt)
    tv <- .isTransversion(ref, alt)
    if (by.sample) {
      isVar <- function(x) {colSums(x, na.rm=TRUE) > 0}
      tisum <- integer(.nSamp(gdsobj))
      tvsum <- integer(.nSamp(gdsobj))
      seqParallel(, gdsobj, FUN = function(f) {
      seqApply(f, "genotype",as.is="none",var.index="relative",
               function(index, x) {
                 tisum <<- tisum + (ti[index] & isVar(x));
                 tvsum <<- tvsum + (tv[index] & isVar(x))
               })},split="by.variant")
      titv <- data.frame(tisum,tvsum)
      if (use.names)
        row.names(titv) <- seqGetData(gdsobj, "sample.id")
      titv
    } else {
      sum(ti) / sum(tv)
    }
  }


sampleQC<-function(gdsfile, sampleid,varlist){

##### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)
print(length(seqGetData(gds, "variant.id")))

#### filter to QCed variants and F9 samples
#varfile<-paste0("/medpop/afib/projects/2020_WGS_WES_ECG/Sequence_QC/TopMed/variantQC/varlist/TOPMed_Freeze9_variantQC_noLCR_Callrate95p_HWE_chr1.tsv")
#vardata<-fread(varfile,header=F,sep="\t",data.table=F)
#phenfile<-"/medpop/afib/projects/2020_WGS_WES_ECG/Phenotypes/grant_tables/data/topmed_ecg_wgs_f9_ids.txt"
#phendata<-fread(phenfile,header=T,data.table=F)
phenlist0<-sampleid
varlist0<-varlist
seqSetFilter(gds, variant.id=varlist0, sample.id=phenlist0)


#### call rate
print("call rate")
samid0 <- seqGetData(gds, "sample.id")
missingrate0<-seqMissing(gds, per.variant=FALSE, .progress=TRUE)
names(missingrate0)<-samid0

##### SNP  Indel ratio
print("SNP  Indel ratio")
varid0 <- seqGetData(gds, "variant.id")
snp0<-isSNV(gds)
names(snp0)<-varid0
snpvar<-names(which(snp0==TRUE))
indelvar<-names(which(snp0==FALSE))

##### SNP
print("SNP")
seqSetFilter(gds, variant.id=snpvar, sample.id=phenlist0)
snphet<-parahetcount(gds)
snphom<-parahomcount(gds,allele="alt",margin="by.sample")
nalt<-snphet$het+snphom$hom*2
snpout<-data.frame(snpalt=nalt,nsnp=snphet$nonmiss)
seqResetFilter(gds)

###### INDEL
print("INDEL")
seqSetFilter(gds, variant.id=indelvar, sample.id=phenlist0)
indelhet<-parahetcount(gds)
indelhom<-parahomcount(gds,allele="alt",margin="by.sample")
nalt<-indelhet$het+indelhom$hom*2
indelout<-data.frame(indelalt=nalt,nindel=indelhet$nonmiss)
seqResetFilter(gds)

##### HET nonrefhom count
print("HET nonrefhom count")
seqSetFilter(gds, variant.id=varlist0, sample.id=phenlist0)
allhet<-snphet$het+indelhet$het
allhom<-snphom$hom+indelhom$hom
nvar<-snphet$nonmiss+indelhet$nonmiss
hethomout<-data.frame(allhet,allhom,nvar)

##### TI TV ratio
print("TI TV ratio")
titvcount<-paratitvcount(gds, by.sample=TRUE,use.names=TRUE)

###### singleton
print("singleton")
ac1count<-countSingletons(gds, use.names=TRUE)

##### total variants included in a chromosome
print("total variants included in a chromosome")
totalvar<-length(varid0)

###### create data set
print("create data set")
out1<-data.frame(sample.id=samid0,missing=missingrate0,snpout,indelout,hethomout,titvcount,ac1count,totalvar)


#### close Info
print("close Info")
seqClose(gds)

#### return result
print("return result")
return(out1)
}


##### call arguments
args <- commandArgs(trailingOnly=T)
gdsfile <- args[1]
nvars<-as.numeric(args[2])
ncpus <- as.numeric(args[3])


#nvars=500
#gdsfile<-"/mnt/project/exome_450k_plink/merged/ukb23156_c10_genotype_variant_QCed_merged.gds"
#ncpus=1
#seqSetFilter(gds, variant.id=varlist0, sample.id=phenlist0)
install.packages(c("doSNOW","doMC"),repos="https://cloud.r-project.org")
BiocManager::install(c('SeqArray',"GENESIS","SeqVarTools"),ask=FALSE,update=FALSE)
##### setup library
suppressMessages(library(SeqVarTools))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(doSNOW))
suppressMessages(library(doMC))
suppressMessages(library(foreach))
cluster = makeCluster(ncpus, type = "SOCK")
registerDoSNOW(cluster)





gds <- seqOpen(gdsfile, allow.duplicate=T)
sampleid <- seqGetData(gds, "sample.id")
varid0 <- seqGetData(gds, "variant.id")
nvar0<-length(varid0)
varidnums<-c(seq(1,nvar0,by=nvars),varid0[nvar0]+1)
totalloops=length(varidnums)-1
seqClose(gds)


basename<-gsub(".gds","",basename(gdsfile))

foreach(ii = 1:totalloops) %dopar% {
#foreach(ii = 1:3) %dopar% {
#
suppressMessages(library(SeqVarTools))
suppressMessages(library(parallel))
suppressMessages(library(data.table))

varlist<-varid0[varidnums[ii]:(varidnums[ii+1]-1)]
out1<-sampleQC(gdsfile, sampleid,varlist)

fwrite(out1,paste0(basename,"_chuck",ii,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")

}

#varidlist_correct<-c(seq(1,tail(varid0, n=1),by=50000),tail(varid0, n=1))


sessionInfo()
quit("no")
