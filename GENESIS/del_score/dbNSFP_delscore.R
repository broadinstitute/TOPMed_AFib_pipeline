####
#### LOF carrier status
files.sources = list.files(path="TOPMed_AFib_pipeline/R/")
files.sources = paste0("TOPMed_AFib_pipeline/R/",files.sources)
sapply(files.sources, source)

##### argument : chromosome info
args=(commandArgs(TRUE))
infile=as.character(args[1])

##### input files
filename<-infile
outfile<-gsub(".gz",".score.tsv",filename)

##### perfrom scoring
result0<-delscore(filename=filename)
write.table(result0,outfile,col.names=T,row.names=F,quote=F,sep="\t")

##### DONE
sessionInfo()
quit("no")
