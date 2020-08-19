#### check wether my wdl script is correct
java -jar /Users/schoi/cromwell/womtool-52.jar validate GENESIS_collapsed_test.wdl

#### what are my input files
java -jar /Users/schoi/cromwell/womtool-52.jar inputs GENESIS_collapsed_test.wdl > GENESIS_collapsed_test.wdl.json






###### example making in R
num=21
gdsfile<-paste0("/Volumes/medpop_afib/data/NHLBI_WGS/TOPMed_phaseI/ncbi/dbGaP-11854/75290/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.8/minDP10.gds/freeze.8.chr",num,".pass_and_fail.gtonly.minDP10.gds")
varfile<-paste0("/Volumes/medpop_afib/schoi/projects/TOPMed/Freeze8/Result/QC/variantQC/varlist/rm_mono/TOPMed_Freeze8_variantQC_noLCR_Callrate95p_HWE_nomono_chr",num,".tsv")
groupfile<-paste0("/Volumes/medpop_afib/data/NHLBI_WGS/TOPMed_phaseI/ncbi/dbGaP-11854/74461/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.8_annotation/var_grouping/hclof_v2/f8_hclof_noflags_chr_",num,"_v2_test.RData")

dat0<-data.frame(num,gdsfile,varfile,groupfile,stringsAsFactors=F)
write.table(dat0,"/Users/schoi/github/TOPMed_AFib_pipeline/GENESIS/collapse/example_inputfile.tsv",col.names=F,row.names=F,quote=F,sep="\t")
