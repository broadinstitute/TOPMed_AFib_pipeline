#### groupfile<-paste0("mnt/f8_hclof_missense_chr_",num,".RData")
#### gdsfile<-paste0("mnt/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10_test.gds")
#### varfile<-paste0("mnt/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,"_test.tsv")
#### phenfile<-"mnt/TOPMed_Freeze8_GQ_AFib_pheno_noMESA_QCed.tsv"
#### nullfile<-"mnt/AF_maleassocatedPCs_kinship.RData"
#### scorefile<-paste0("mnt/TOPMed_Freeze8_missense_chr",num,".annot.vcf.score.tsv")
#### scorename<-c("Dprop")
#### scutoff<-0.9
#### acutoff<-0.001
#### outfile<-paste0("mnt/TOPMed_freeze8_AF_hclof_missense_noweight_",paste0(scorename,collapse="_"),"_30delscore_noweighted_cutoff",scutoff,"_chr",num,".tsv")
#### stat="Burden"

####
#### LOF carrier status
files.sources = list.files(path="TOPMed_AFib_pipeline/R/")
files.sources = paste0("TOPMed_AFib_pipeline/R/",files.sources)
sapply(files.sources, source)

##### argument : chromosome info
args=(commandArgs(TRUE))
infile=as.character(args[1])

num=as.character(args[1])
gdsfile=as.character(args[2])
varfile=as.character(args[3])
groupfile=as.character(args[4])
phenfile=as.character(args[5])
nullfile=as.character(args[6])
scorefile=as.character(args[7])
scorename=as.character(args[8])
scutoff=as.numeric(args[9])
acutoff=as.numeric(args[10])
stat=as.character(args[11])
outfile=as.character(args[12])


results<-del.cutoff.Burden(num,gdsfile,groupfile,scorefile,scorename,scutoff,acutoff,phenfile,nullfile,stat,outfile)


##### DONE
sessionInfo()
quit("no")
