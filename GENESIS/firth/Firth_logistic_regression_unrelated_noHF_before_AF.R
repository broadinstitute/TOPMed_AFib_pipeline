####
#### LOF carrier status
files.sources = list.files(path="TOPMed_AFib_pipeline/R/")
files.sources = paste0("TOPMed_AFib_pipeline/R/",files.sources)
sapply(files.sources, source)

##### argument : chromosome info
args=(commandArgs(TRUE))
phenfile=as.character(args[1])
num=as.character(args[2])
gdsfile=as.character(args[3])
genefile=as.character(args[4])
resultfile=as.character(args[5])
unrelcol=as.character(args[6])
noHF_col=as.character(args[7])
outfile=as.character(args[8])

####
#### loop by the geneid
genelist<-fread(genefile,header=F,sep="\t",data.table=F)
chrgenes<-subset(genelist,chr==num)
geneids<-chrgenes$geneid

result<-list()
for (gnum in 1:length(geneids)){
geneid<-geneids[gnum]

#### perfrom the firth
mod1<-unrelated_firth(phenfile=phenfile,num=num,gdsfile=gdsfile,geneid=geneid,resultfile=resultfile,unrelcol=unrelcol,noHF_col=noHF_col)
result[[geneid]]<-mod1
}
save(result,filename=outfile)
