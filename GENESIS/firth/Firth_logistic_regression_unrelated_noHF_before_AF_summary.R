#####
#####
files.sources = list.files(path="TOPMed_AFib_pipeline/R/")
files.sources = paste0("TOPMed_AFib_pipeline/R/",files.sources)
sapply(files.sources, source)


##### argument : chromosome info
args=(commandArgs(TRUE))
print(args)
outfiles<-unlist(strsplit(args[1], ","))

#######
####### read files
splits0<-gsub("Chr","",unlist(strsplit(basename(outfiles),"_",fixed=T)))
chrs<-splits0[seq(1,length(splits0),by=4)]
print(outfiles)
print(chrs)

combresult<-list()
for (ii in 1:length(chrs)){

num<-chrs[ii]
outfile<-outfiles[ii]
result0<-get(load(outfile))
gids<-names(result0)
for (jj in 1:length(gids)){

combresult[[gids[jj]]]<-result0[[gids[jj]]]
}
}
save(combresult,file="summary.RData")



sessionInfo()
quit("no")
