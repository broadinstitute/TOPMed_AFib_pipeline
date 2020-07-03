#######
#######
cd /home/jupyter-user/notebooks/BioData_catalyst_afib/edit/TOPMed_AFib_pipeline
git pull https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

######## copy dbnsfp file from gc
gsutil cp gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/v1/chr*_dbNSFP_scatter_piece.gathered.sorted.tsv.gz /home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/



#########
source("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/TOPMed_AFib_pipeline/R/dbnsfpfilter.R")

for ( num in c(1:22,"X")){

####### dbnsfp dataset
filename<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/chr",num,"_dbNSFP_scatter_piece.gathered.sorted.tsv.gz")
normvar<-c("chr","pos","ref","alt","Ensembl_geneid","Ensembl_transcriptid","MetaLR_pred","MetaSVM_pred","MutationAssessor_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred",
"PROVEAN_pred","M_CAP_pred","FATHMM_pred","SIFT4G_pred","SIFT_pred","PrimateAI_pred","DEOGEN2_pred")
scorevar<-c("MutPred_rankscore","REVEL_rankscore","VEST4_rankscore")
outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/s0/chr",num,"_dbNSFP_scatter_piece.gathered.sorted.pred.tsv")

###### write table
res<-dbnsfpscore(filename=filename,normvar=normvar,scorevar=scorevar,phredvar=NULL)
}

#########
#########
cpcmd<-"gzip /home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/s0/chr*_dbNSFP_scatter_piece.gathered.sorted.pred.tsv"
system(cpcmd,intern=T)
cpcmd<-"gsutil cp /home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/s0/chr*_dbNSFP_scatter_piece.gathered.sorted.pred.tsv.gz gs://fc-e6874e4a-c83c-4c52-b370-1da3c49153fe/annot/s0/"
system(cpcmd,intern=T)



##########
##########
for (num in c(1:22,"X")){

outfile<-paste0("/home/jupyter-user/notebooks/BioData_catalyst_afib/edit/annot/wgsa/s0/chr",num,"_dbNSFP_scatter_piece.gathered.sorted.pred.tsv")

dat1<-fread(cmd=paste0("gunzip -c ",outfile),header=T,data.table=F,sep="\t")
dat1$varid<-paste(dat1$chr,dat1$pos,dat1$ref,dat1$alt,dat1$Ensembl_geneid,sep=":")
dat2<-subset(dat1,nmiss>=5)
dat3<-dat2[order(dat2$varid,-dat2$dprop),]
dat4<-dat3[!duplicated(dat3$varid),]
}
