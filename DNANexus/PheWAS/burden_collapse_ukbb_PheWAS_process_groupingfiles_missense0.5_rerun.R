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
  old_file <- paste0('/mnt/project/sjj/projects/phewas/v1/results/association/', phecode, '_results_chr', chr, '_maf0.00001.RData')
  new_file <- paste0('/mnt/project/sjj/projects/phewas/v1/data/grouping_files/ukbb_phewas_v1_groupingfile_c', chr, '_hclofnoflagmissense0.5_missense0.2_popmax0.00001.RData')
  
  old_results <- get(load(old_file))
  group <- get(load(new_file))
  group <- group[which(grepl("missense0.2", group$group_id)), ] 
  group <- group[group$Dprop>=0.5 & group$Dtools>=7, ]
  group$group_id <- gsub("missense0.2", "missense0.5", group$group_id)
  group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)
  
  new_group <- NULL
  #new_res1 <- old_results
  genes <- unique(gsub("_missense0.5.*", "", rownames(old_results$results[old_results$results$n.alt>0, ])))
  
  for(gene in genes){
      inter_group <- group[which(grepl(gene, group$group_id)), ]
      new_vars <- inter_group$varid
      old_vars <- old_results$variantInfo[[paste0(gene, '_missense0.5_gnomAD_POPMAX0.00001')]]
      old_vars$varid <- paste0(old_vars$chr, ":", old_vars$pos, ":", old_vars$ref, ":", old_vars$alt)
      old_vars <- old_vars$varid
  
      #check
      if(!all(old_vars %in% new_vars)){
          new_group <- rbind(new_group, inter_group)
      }
  }
  
  cat(length(unique(new_group$group_id)), " out of ", length(genes), " equals ", length(unique(new_group$group_id))/length(genes)*100 , "%.\n")
  
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
  lowmem_vars <- group_tally[which(group_tally[,2]<100), 1]
  #length(lowmem_vars)
  save(new_group[new_group$group_id %in% lowmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_missense0.5_popmax0.00001_lowmem.RData'))
  
  ### highmem
  highmem_vars <- group_tally[group_tally[,2]>=100 & group_tally[,2]<200, 1]
  #length(highmem_vars)
  save(new_group[new_group$group_id %in% highmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_missense0.5_popmax0.00001_highmem.RData'))
  
  ### highhighmem
  highhighmem_vars <- group_tally[group_tally[,2]>=200 & group_tally[,2]<300, 1]
  #length(highhighmem_vars)
  save(new_group[new_group$group_id %in% highhighmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_missense0.5_popmax0.00001_highhighmem.RData'))
  
  ### veryhighmem
  veryhighmem_vars <- group_tally[group_tally[,2]>=300, 1]
  #length(veryhighmem_vars)
  save(new_group[new_group$group_id %in% veryhighmem_vars, ], 
      save=paste0('ukbb_phewas_v1_groupingfile_c', chr, '_missense0.5_popmax0.00001_veryhighmem.RData'))
  
}

sessionInfo()
quit("no")
