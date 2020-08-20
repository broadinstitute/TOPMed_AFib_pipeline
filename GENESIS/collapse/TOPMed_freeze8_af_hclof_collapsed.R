#####
#####
##### R.utils::sourceDirectory("TOPMed_AFib_pipeline/R/", modifiedOnly=TRUE);
files.sources = list.files(path="TOPMed_AFib_pipeline/R/")
files.sources = paste0("TOPMed_AFib_pipeline/R/",files.sources)
sapply(files.sources, source)


##### argument : chromosome info
args=(commandArgs(TRUE))
num=as.character(args[1])
gdsfile=as.character(args[2])
varfile=as.character(args[3])
groupfile=as.character(args[4])
phenfile=as.character(args[5])
nullfile=as.character(args[6])
stat=as.character(args[7])
cutoff=as.numeric(args[8])
outfile=as.character(args[9])

##### function start
res0<-aggreated_test(num=num,gdsfile=gdsfile,varfile=varfile,groupfile=groupfile,phenfile=phenfile,nullfile=nullfile,cutoff=cutoff,stat=stat,outfile=outfile)


sessionInfo()
quit("no")
