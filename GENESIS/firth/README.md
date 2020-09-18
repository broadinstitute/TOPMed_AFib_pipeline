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
for (num in c(1:22,"X")){

gdsfile<-paste0("/Users/schoi/github/testdata/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10_test.gds")
varfile<-paste0("/Users/schoi/github/testdata/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,"_test.tsv")
groupfile<-paste0("/Users/schoi/github/testdata/f8_hclof_missense_chr_",num,".RData")
dbnsfpfile<-"/Users/schoi/github/testdata/TOPMed_Freeze8_missense_chr22.annot.vcf.gz"
scorefile<-"TOPMed_Freeze8_missense_chr22.annot.vcf.score.tsv"
dat1<-data.frame(num,gdsfile,varfile,groupfile,dbnsfpfile,scorefile,stringsAsFactors=F)
dat0<-rbind(dat0,dat1)
}
write.table(dat0,"/Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/del_score/inputfile_chrall.tsv",col.names=F,row.names=F,quote=F,sep="\t")




### add drsadress
library(jsonlite)
setwd("/medpop/afib/schoi/projects/TOPMed/Freeze8/Result/terra")
dat0<-NULL
jsonread<-read_json("dsr_gds_file_manifest.json",TRUE)
rownames(jsonread)<-jsonread[,"file_name"]
jsonread2<-jsonread[paste0("freeze.8.chr",c(1:22,"X"),".pass_and_fail.gtonly.minDP10.gds"),]
nums<-c(1:22,"X")
for (ii in c(1:23)){

num<-nums[ii]
gdsfile<-paste0("drs://",jsonread2[ii,"object_id"])
varfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/gds/varlist/rm_mono/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
groupfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/var_grouping/hclof_missense_clinvar_v2/f8_hclof_noflags_missense_clinvarPLP_chr_",num,".RData")
dbnsfpfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/dbnsfp/delscore30/TOPMed_Freeze8_missense_chr",num,".annot.vcf.gz")
scorefile<-paste0("TOPMed_Freeze8_missense_chr",num,".annot.vcf.score.tsv")

dat1<-data.frame(num,gdsfile,varfile,groupfile,dbnsfpfile,scorefile,stringsAsFactors=F)
dat0<-rbind(dat0,dat1)
}
write.table(dat0,"inputfile_GEN3_hclofnoflag_missense0.9_clinvarPLP_chrall.tsv",col.names=F,row.names=F,quote=F,sep="\t")


cd /medpop/afib/schoi/projects/TOPMed/Freeze8/Result/terra
gsutil cp inputfile_GEN3_hclofnoflag_missense0.9_clinvarPLP_chrall.tsv gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/WDL_inputs/

gsutil cp gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/var_grouping/hclof_missense_clinvar_v2/f8_hclof_noflags_missense_clinvarPLP_chr_1.RData .

#### significant genes in EUR

### add drsadress
library(jsonlite)
setwd("/medpop/afib/schoi/projects/TOPMed/Freeze8/Result/terra")
dat0<-NULL
jsonread<-read_json("dsr_gds_file_manifest.json",TRUE)
rownames(jsonread)<-jsonread[,"file_name"]
jsonread2<-jsonread[paste0("freeze.8.chr",c(1:22,"X"),".pass_and_fail.gtonly.minDP10.gds"),]
nums<-c(1:22,"X")
for (ii in c(1,2,6,8,10,11,17)){

num<-nums[ii]
gdsfile<-paste0("drs://",jsonread2[ii,"object_id"])
varfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/gds/varlist/rm_mono/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
groupfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/var_grouping/hclof_missense_v2/f8_hclof_noflags_missense_suggestive_chr_",num,"_v2.RData")
dbnsfpfile<-paste0("gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/dbnsfp/delscore30/TOPMed_Freeze8_missense_chr",num,".annot.vcf.gz")
scorefile<-paste0("TOPMed_Freeze8_missense_chr",num,".annot.vcf.score.tsv")

dat1<-data.frame(num,gdsfile,varfile,groupfile,dbnsfpfile,scorefile,stringsAsFactors=F)
dat0<-rbind(dat0,dat1)
}
write.table(dat0,"inputfile_GEN3_hclofnoflag_missense0.9_suggestive.tsv",col.names=F,row.names=F,quote=F,sep="\t")

cd /medpop/afib/schoi/projects/TOPMed/Freeze8/Result/terra
gsutil cp inputfile_GEN3_hclofnoflag_missense0.9_suggestive.tsv gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/WDL_inputs/





#### significant genes in no HF before AF
