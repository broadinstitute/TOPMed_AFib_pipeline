##### manoutfile<-paste0("/medpop/afib/schoi/projects/TOPMed/Freeze8/Result/association/rare/hclof/TOPMed_freeze8_af_hclof_burden_manhattan.png")
##### qqoutfile<-paste0("/medpop/afib/schoi/projects/TOPMed/Freeze8/Result/association/rare/hclof/TOPMed_freeze8_af_hclof_burden_qq.png")
gencodgene<-"TOPMed_AFib_pipeline/GENESIS/data/gencode_v28_genes.gz"

#####
#####
##### R.utils::sourceDirectory("TOPMed_AFib_pipeline/R/", modifiedOnly=TRUE);
files.sources = list.files(path="TOPMed_AFib_pipeline/R/")
files.sources = paste0("TOPMed_AFib_pipeline/R/",files.sources)
sapply(files.sources, source)
library(qqman)


##### argument : chromosome info
args=(commandArgs(TRUE))
print(args)
maccutoff<-as.numeric(args[1])
outfiles<-unlist(strsplit(input_args[2], ","))
##### gencodgene=as.character(args[1])
##### manoutfile=as.character(args[2])
##### qqoutfile=as.character(args[3])


#######
####### read files
chrs<-gsub("Chr","",gsub("_collapsed_results.RData","",outfiles))
print(outfiles)
print(chrs)
sumres0<-summarydata(files=outfiles,chrs=chrs,thre_cMAC=maccutoff)
result0<-sumres0$generesult
result0$gene<-rownames(result0)

######## protein coding genes
######## merge with protein coding genes
genes<-fread(cmd=paste0("zcat ",gencodgene),header=T,data.table=F,sep="\t")
pcgenes<-subset(genes,gene_type=="protein_coding" & chr!="Y")

########
######## merge with protein coding genes
result1<-merge(result0,pcgenes,by.x="gene",by.y="id")
p_thre<-0.05/(nrow(result1))

########
########
save(result1,file="summary.RData")

########
######## make a manhattan plit
png("manhattan_plot.png", width=1500,height=500,type="cairo",res=100)
par(mar = c(5.1, 6.1, 5.1, 2.1))
man(result1, chr="chr", bp="mpos", snp="gene", p="Score.pval",
         chrlabs=c(1:22, "X"),
         main="",
         col=c("dodgerblue4", "firebrick4"),
         genomewideline= -log10(p_thre),
         suggestiveline = -log10(0.1/(nrow(result1))),
         cex.lab=1.5)
axis(side=2,at=seq(1,16,by=2),labels=seq(1,16,by=2), las = 1, pch = 20,cex.axis=0.8,tck=-0.03)
dev.off()

######
###### make a qqplot
lambda <- calculateLambda((result0$Score.Stat)^2, df=1)
png("qqplot.png", width=1000,height=1000,type="cairo",res=100)
qq(result0$Score.pval)
legend("topleft",bty = "n",legend=paste0("Lambda = ",round(lambda,digits=2)))
dev.off()


sessionInfo()
quit("no")
