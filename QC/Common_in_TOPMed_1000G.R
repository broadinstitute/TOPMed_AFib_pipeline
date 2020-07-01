#######
#######
/tmp/plink2

gsutil cp gs://fc-306d0fc4-2f1d-4ea1-ae29-ea5b8fe0cb22/QC/High_Quality_1000G/ALL.chr10_GRCh38.genotypes.20170504.cleaned.* ./


grep -v ^# ALL.chr10_GRCh38.genotypes.20170504.cleaned.pvar | awk -F"\t" '{print $1"\t"$2"\t"$1":"$2":"$4":"$5"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | head -n 1 > test
for ii in {1..22}
do
grep ^# ALL.chr${ii}_GRCh38.genotypes.20170504.cleaned.pvar > newpvar/ALL.chr${ii}_GRCh38.genotypes.20170504.cleaned.pvar
grep -v ^# ALL.chr${ii}_GRCh38.genotypes.20170504.cleaned.pvar | awk -F"\t" '{print $1"\t"$2"\t"$1":"$2":"$4":"$5"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' >> newpvar/ALL.chr${ii}_GRCh38.genotypes.20170504.cleaned.pvar
done

cd /medpop/afib/data/1000Genomes/phase3_v5a/Build38/cleaned/clean0/newpvar
gsutil cp *.* gs://fc-306d0fc4-2f1d-4ea1-ae29-ea5b8fe0cb22/QC/High_Quality_1000G/


######
gsutil cp gs://fc-306d0fc4-2f1d-4ea1-ae29-ea5b8fe0cb22/QC/High_Quality_1000G/ALL.chr10_GRCh38.genotypes.20170504.cleaned.* ./

####### R

##### export information
for(num in c(1:22)){

topmedfile<-paste0("TOPMed_Freeze8_chr",num,"_Highquality_variants.tsv")
##### copy result to the bucket
com0<-paste0("gsutil cp  ",bucket,"/QC/High_Quality_TOPMed/",topmedfile," ./")
system(com0)
##### high qulity topmed
hqtopmed<-fread(topmedfile,header=T,data.table=F,sep="\t")



##### high qulity 1000G
tgfile<-paste0("ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned.pvarID")
##### copy result to the bucket
com0<-paste0("gsutil cp  ",bucket,"/QC/High_Quality_1000G/",tgfile," ./")
system(com0)
###### read High quality 1000G
hq1000g<-fread(tgfile,header=F,data.table=F,sep="\t")

###### common in 1000G
hq_1000g<-hq1000g[hq1000g$V1 %in% hqtopmed$varid,]
hq_topmed<-hqtopmed[hqtopmed$varid %in% hq1000g$V1,]


###### copy TOPMed file to cloud
topmedfile1<-gsub(".tsv","_in_comm.tsv",topmedfile)
write.table(hq_topmed,topmedfile1,col.names=T,row.names=F,quote=F,sep="\t")
##### copy result to the bucket
com0<-paste0("gsutil cp  ",topmedfile1," ",bucket,"/QC/High_Quality_TOPMed/Incommon/")
system(com0)


###### copy 1000G file to cloud
tgfile1<-gsub(".pvarID","_in_comm.tsv",tgfile)
write.table(hq_1000g,tgfile1,col.names=F,row.names=F,quote=F,sep="\t")
##### copy result to the bucket
com0<-paste0("gsutil cp  ",tgfile1," ",bucket,"/QC/High_Quality_1000G/Incommon/")
system(com0)



print(paste0(("chromosome",num,"is done"))
}



for ( num in 1:22){

##### high qulity 1000G
tgfile<-paste0("ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned.*")
##### copy result to the bucket
com0<-paste0("gsutil cp  ",bucket,"/QC/High_Quality_1000G/",tgfile," ./")
system(com0)
###### read High quality 1000G
tgfile0<-gsub(".*","",tgfile,fixed=T)
tgfile1<-gsub(".*","_in_comm.tsv",tgfile,fixed=T)
tgfile2<-gsub(".*","_in_comm",tgfile,fixed=T)

####### select variants
com1<-paste0("/tmp/plink2/plink2 --pfile ",tgfile0," --extract ",tgfile1," --make-bed --out ",tgfile2)
system(com1,intern=T)

########
######## copy files to bucket
tgfile3<-gsub("_in_comm","_in_comm.*",tgfile2,fixed=T)

com2<-paste0("gsutil cp  ",tgfile3," ",bucket,"/QC/High_Quality_1000G/Incommon/")
system(com2,intern=T)
print(paste0("chromosome",num," is done"))

}

######
sessionInfo()
quite("no")


#########
######### perform link
##### high qulity 1000G
##### copy result to the bucket

######
###### Copy superpopulation file
bucket<- Sys.getenv('WORKSPACE_BUCKET')
tgfile<-paste0("subpop1000G.*")
copy1000Gpop<-paste0("gsutil cp  ",bucket,"/QC/High_Quality_1000G/superpop/",tgfile," ./")
system(copy1000Gpop,intern=T)

######
######  Copy plink formatted file to server
for ( num in 1:22){

#####
tgfile<-paste0("ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.*")
copy1000Gbed<-paste0("gsutil cp  ",bucket,"/QC/High_Quality_1000G/Incommon/",tgfile," ./")
system(copy1000Gbed,intern=T)

prEUR<-paste0("/tmp/plink2/plink2 --bfile ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm --keep subpop1000G.EUR.txt --indep-pairwise 50 5 0.2 --out ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR")
prEAS<-paste0("/tmp/plink2/plink2 --bfile ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm --keep subpop1000G.EAS.txt --extract ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.prune.in --indep-pairwise 50 5 0.2 --out ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS")
prSAS<-paste0("/tmp/plink2/plink2 --bfile ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm --keep subpop1000G.SAS.txt --extract ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.prune.in --indep-pairwise 50 5 0.2 --out ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS")
prAMR<-paste0("/tmp/plink2/plink2 --bfile ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm --keep subpop1000G.AMR.txt --extract ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.prune.in --indep-pairwise 50 5 0.2 --out ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.AMR")
prAFR<-paste0("/tmp/plink2/plink2 --bfile ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm --keep subpop1000G.AFR.txt --extract ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.AMR.prune.in --indep-pairwise 50 5 0.2 --out ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.AMR.AFR")
prAll<-paste0("/tmp/plink2/plink2 --bfile ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm --extract ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.AMR.AFR.prune.in --make-bed --out ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.pruned")

system(prEUR,intern=T)
system(prEAS,intern=T)
system(prSAS,intern=T)
system(prAMR,intern=T)
system(prAFR,intern=T)
system(prAll,intern=T)

}
system(paste0("wc -l  ALL.chr*_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.AMR.AFR.prune.in"),intern=T)


#########
######### read gds file and select the variants
library(data.table)
library(SeqArray)

for (num in 1:22){

###### copy common variants
###### copy TOPMed file to cloud
topmedfile<-paste0("TOPMed_Freeze8_chr",num,"_Highquality_variants_in_comm.tsv")
copytopmedfile<-paste0("gsutil cp ",bucket,"/QC/High_Quality_TOPMed/Incommon/",topmedfile," ./")
system(copytopmedfile)
tmcom<-fread(topmedfile,header=T,sep="\t")

######
###### 1000G pruned variant list
tgfile<-paste0("ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.EUR.EAS.SAS.AMR.AFR.prune.in")
tgcom<-fread(tgfile,header=F,sep="\t")

#######
####### pruned list
tmcompr<-tmcom[tmcom$varid %in% tgcom$V1,]
varlist<-tmcompr$gdsid

#######
####### copy gds file
gdsfile<-paste0("freeze.8.chr",num,".pass_only.phased.gds")
copycommand<-paste0("gsutil -m cp ",bucket,"/genotype/freeze8_gds/",gdsfile," ./")
system(copycommand)

####### read gds file
gds <- seqOpen(gdsfile, allow.duplicate=T)

####### filered gds file
seqSetFilter(gds, variant.id=varlist)

####### export selected variants to vcf format
vcffile<-gsub(".gds",".pruned.vcf.gz",gdsfile)
seqGDS2VCF(gds, vcffile, info.var=character(),fmt.var=character())
seqClose(gds)
}

##########
########## combine TOPMed and 1000G plink files by chromosome

for ( num in 1:22){
vcffile<-paste0("freeze.8.chr",num,".pass_only.phased.pruned.vcf.gz")
bfile<-paste0("ALL.chr",num,"_GRCh38.genotypes.20170504.cleaned_in_comm.pruned")
outfile<-paste0("freeze.8.chr",num,".pass_only.phased.pruned.1000G.combined")
vcftoplink<-gsub(".vcf.gz","",vcffile)
### step1
chrplink<-paste0("/tmp/plink/plink --vcf ",vcffile," --keep-allele-order --make-bed --out ",vcftoplink)
system(chrplink)

##### fix the variant ID
fixid<-paste0("awk '{print $1\"\t\"$1\":\"$4\":\"$6\":\"$5\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}' ",vcftoplink,".bim > ",vcftoplink,".new")
system(fixid)
mvfile<-paste0("mv ",vcftoplink,".new ",vcftoplink,".bim")
system(mvfile)

#####
chrcombine<-paste0("/tmp/plink/plink --bfile ",vcftoplink," --bmerge ",bfile," --keep-allele-order --out ",outfile)
system(chrcombine)
}

##########
########## combine all chromosomes

beds<-paste0("freeze.8.chr",2:22,".pass_only.phased.pruned.1000G.combined.bed")
bims<-paste0("freeze.8.chr",2:22,".pass_only.phased.pruned.1000G.combined.bim")
fams<-paste0("freeze.8.chr",2:22,".pass_only.phased.pruned.1000G.combined.fam")
mlist<-data.frame(beds,bims,fams)
write.table(mlist,"plink_merge_list.txt",col.names=F,row.names=F,quote=F,sep="\t")

########## combine
chr1file<-"freeze.8.chr1.pass_only.phased.pruned.1000G.combined"
outfile<-"freeze.8.chrall.pass_only.phased.pruned.1000G.combined"
combinechr<-paste0("/tmp/plink/plink --bfile ",chr1file," --merge-list plink_merge_list.txt --keep-allele-order --make-bed --out ",outfile)
system(combinechr)

##### copy result to the bucket
outfile<-"freeze.8.chrall.pass_only.phased.pruned.1000G.combined.*"
copyprunedfile<-paste0("gsutil cp ",outfile," ",bucket,"/QC/High_Quality_Combined/")
system(copyprunedfile,intern=T)
