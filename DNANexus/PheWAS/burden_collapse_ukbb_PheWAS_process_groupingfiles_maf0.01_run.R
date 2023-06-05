#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
print(args)
gdsfile=as.character(args[1])
groupfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
outfile=as.character(args[6])
AF.max=as.numeric(args[7])
MAC.max=as.numeric(args[8])
score.method=as.character(args[9])
key.file=as.character(args[10])
key.file.phecodecol=as.character(args[11])
phenum=as.numeric(args[12])
cutoff1=as.numeric(args[13])
cutoff2=as.numeric(args[14])
cutoff3=as.numeric(args[15])

.libPaths(c("rpackages4_1_3",.libPaths()))

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("TOPMed_AFib_pipeline/DNANexus/kernell_variance_component_modfied.R")

key <- fread(key.file, stringsAsFactors=F, data.table=F)
phecodes <- unique(key[,key.file.phecodecol])
phecode <- phecodes[phenum]
nullfile <- gsub("PHENUM", phecode, nullfile)
outfile <- gsub("PHENUM", phecode, outfile)

chr <- gsub("_genotype_variant_sample_QCed.gds", "", gsub("tmp/ukb23156_c", "", gdsfile))
memlevel <- gsub(".RData", "", gsub(".*_", "", groupfile))
outfile_saved <- paste0("/mnt/project/sjj/projects/phewas/v1/results/association/round2/", memlevel, "/chr", chr, "/", outfile)

if(!file.exists(outfile_saved)){
  new_file <- paste0('/mnt/project/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c', chr, '_hclofnoflagmissense0.8missense0.5_popmax0.01.RData')
  
  group <- get(load(new_file))
  group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
  
  new_group <- group
  
  ### Split by mem
  groups <- unique(new_group$group_id)
  group_tally <- NULL
  for(group in groups){
      N_vars <- nrow(new_group[new_group$group_id==group, ])
      group_tally <- rbind(group_tally, c(group, N_vars))
  }
  group_tally <- as.data.frame(group_tally)
  group_tally[,2] <- as.numeric(group_tally[,2])
  
  ### lowmem
  lowmem_vars <- group_tally[which(group_tally[,2]<cutoff1), 1]
  #length(lowmem_vars)
  save(new_group[new_group$group_id %in% lowmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_lowmem.RData'))
  
  ### highmem
  highmem_vars <- group_tally[group_tally[,2]>=cutoff1 & group_tally[,2]<cutoff2, 1]
  #length(highmem_vars)
  save(new_group[new_group$group_id %in% highmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highmem.RData'))
  
  ### highhighmem
  highhighmem_vars <- group_tally[group_tally[,2]>=cutoff2 & group_tally[,2]<cutoff3, 1]
  #length(highhighmem_vars)
  save(new_group[new_group$group_id %in% highhighmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_hclofnoflag_hclofnoflagmissense0.8_hclofnoflagmissense0.5_popmax0.01_highhighmem.RData'))
  
  ### veryhighmem
  veryhighmem_vars <- group_tally[group_tally[,2]>=cutoff3, 1]
  #length(veryhighmem_vars)
  save(new_group[new_group$group_id %in% veryhighmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_hclofnoflagmissense0.8missense0.5_popmax0.01_veryhighmem.RData'))
  
}

sessionInfo()
quit("no")
