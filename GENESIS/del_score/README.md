#### check wether my wdl script is correct
java -jar /Users/schoi/cromwell/womtool-52.jar validate GENESIS_delsocre_burden_test.wdl

#### what are my input files
java -jar /Users/schoi/cromwell/womtool-52.jar inputs GENESIS_delsocre_burden_test.wdl > GENESIS_delsocre_burden_test.wdl.json

#### excute!!!!
java -jar /Users/schoi/cromwell/cromwell-52.jar run /Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/del_score/GENESIS_delscore_burden_test.wdl --inputs /Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/del_score/GENESIS_delscore_burden_test.wdl_example.json





###### example making in R
dat0<-NULL
for (num in c(22)){

gdsfile<-paste0("/Users/schoi/github/testdata/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10_test.gds")
varfile<-paste0("/Users/schoi/github/testdata/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,"_test.tsv")
groupfile<-paste0("/Users/schoi/github/testdata/f8_hclof_missense_chr_",num,".RData")
dbnsfpfile<-"/Users/schoi/github/testdata/TOPMed_Freeze8_missense_chr22.annot.vcf.gz"
scorefile<-"TOPMed_Freeze8_missense_chr22.annot.vcf.score.tsv"
dat1<-data.frame(num,gdsfile,varfile,groupfile,dbnsfpfile,scorefile,stringsAsFactors=F)
dat0<-rbind(dat0,dat1)
}
write.table(dat0,"/Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/del_score/inputfile_chr22.tsv",col.names=F,row.names=F,quote=F,sep="\t")




dat0<-NULL
for (num in c(21:22)){

gdsfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/gds/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/gds/varlist/rm_mono/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
groupfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/var_grouping/hclof_v2/f8_hclof_noflags_chr_",num,"_v2.RData")
dbnsfpfile<-
scorefile<-
dat1<-data.frame(num,gdsfile,varfile,groupfile,stringsAsFactors=F)
dat0<-rbind(dat0,dat1)
}
write.table(dat0,"/Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/collapse/inputfile_chr21_22.tsv",col.names=F,row.names=F,quote=F,sep="\t")

####
dat0<-NULL
for (num in c(1:22,"X")){

gdsfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/gds/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/gds/varlist/rm_mono/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
groupfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/var_grouping/hclof_v2/f8_hclof_noflags_chr_",num,"_v2.RData")

dat1<-data.frame(num,gdsfile,varfile,groupfile,stringsAsFactors=F)
dat0<-rbind(dat0,dat1)
}
write.table(dat0,"/Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/collapse/inputfile_chrall.tsv",col.names=F,row.names=F,quote=F,sep="\t")



/medpop/afib/schoi/software/cromwell
