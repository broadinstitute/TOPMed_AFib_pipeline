#for(code_num in c(2:length(codes))){
  #code <- codes[code_num]
args <- commandArgs(trailingOnly = TRUE)
print(args)
chunk_num <- as.numeric(args[1])
n_chunks <- as.numeric(args[2])
rerun <- as.logical(args[3])
if(is.na(rerun) | is.null(rerun)){rerun <- FALSE}

if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}

files1 <- list.files()
files1 <- files1[which(grepl("UKB_sumstats", files1))]
files1 <- files1[which(grepl("formeta", files1))]
length(files1)

files2 <- list.files()
files2 <- files2[which(grepl("AoUv8_sumstats", files2))]
files2 <- files2[which(grepl("formeta", files2))]
length(files2)

files3 <- list.files()
files3 <- files3[which(grepl("MGB_sumstats", files3))]
files3 <- files3[which(grepl("formeta", files3))]
length(files3)

files <- c(files1, files2, files3)
codes <- gsub(".*_sumstats_", "", files)
codes <- gsub("_bothsexes.*", "", codes)
codes <- gsub("_femaleonly.*", "", codes)
codes <- gsub("_maleonly.*", "", codes)
codes <- unique(codes)
message("total codes: ", length(codes))
message("example code: ", codes[1])
#head(codes,n=30)
if(rerun){
  rerun_codes <- rbind(data.table::fread('phecodes_need_removed.tsv', stringsAsFactors=F, data.table=F, header=F),
                       data.table::fread('phecodes_need_rerun.tsv', stringsAsFactors=F, data.table=F, header=F)
                 )
  codes <- codes[which(codes%in%rerun_codes[,1])]
  message("rerun codes: ", length(codes))
  message("example rerun code: ", codes[1])
}

chunks <- split(c(1:length(codes)), cut(seq_along(codes), n_chunks, labels = FALSE))
chunk <- chunks[[chunk_num]]
print(codes[chunk])

for(code_num in chunk){
  code <- codes[code_num]

  message("  ")
  message("  ")
  message("  ")
  message("Busy with ", code, " which is number", code_num)
  message("  ")
  message("  ")
  message("  ")
  
  sex <- 'both_sexes'
  
  # find ukb file
  if(file.exists(paste0('UKB_sumstats_', code, "_bothsexes_mincarriers10_formeta.tsv.gz"))){
    ukb_file <- paste0('UKB_sumstats_', code, "_bothsexes_mincarriers10_formeta.tsv.gz")
  }else if(file.exists(paste0('UKB_sumstats_', code, "_femaleonly_mincarriers10_formeta.tsv.gz"))){
    ukb_file <- paste0('UKB_sumstats_', code, "_femaleonly_mincarriers10_formeta.tsv.gz")
    sex <- 'femaleonly'
  }else if(file.exists(paste0('UKB_sumstats_', code, "_maleonly_mincarriers10_formeta.tsv.gz"))){
    ukb_file <- paste0('UKB_sumstats_', code, "_maleonly_mincarriers10_formeta.tsv.gz")
    sex <- 'maleonly'
  }else{
    ukb_file <- NA
  }
  message("UKB file: ", ukb_file)
  
  # find aou file
  if(file.exists(paste0('AoUv8_sumstats_', code, "_bothsexes_mincarriers20_formeta.tsv.gz"))){
    aou_file <- paste0('AoUv8_sumstats_', code, "_bothsexes_mincarriers20_formeta.tsv.gz")
  }else if(file.exists(paste0('AoUv8_sumstats_', code, "_femaleonly_mincarriers20_formeta.tsv.gz"))){
    aou_file <- paste0('AoUv8_sumstats_', code, "_femaleonly_mincarriers20_formeta.tsv.gz")
    sex <- 'femaleonly'
  }else if(file.exists(paste0('AoUv8_sumstats_', code, "_maleonly_mincarriers20_formeta.tsv.gz"))){
    aou_file <- paste0('AoUv8_sumstats_', code, "_maleonly_mincarriers20_formeta.tsv.gz")
    sex <- 'maleonly'
  }else{
    aou_file <- NA
  }  
  message("AoUv8 file: ", aou_file)
  
  # find mgb file
  if(file.exists(paste0('MGB_sumstats_', code, "_bothsexes_mincarriers10_formeta.tsv.gz"))){
    mgb_file <- paste0('MGB_sumstats_', code, "_bothsexes_mincarriers10_formeta.tsv.gz")
  }else if(file.exists(paste0('MGB_sumstats_', code, "_femaleonly_mincarriers10_formeta.tsv.gz"))){
    mgb_file <- paste0('MGB_sumstats_', code, "_femaleonly_mincarriers10_formeta.tsv.gz")
    sex <- 'femaleonly'
  }else if(file.exists(paste0('MGB_sumstats_', code, "_maleonly_mincarriers10_formeta.tsv.gz"))){
    mgb_file <- paste0('MGB_sumstats_', code, "_maleonly_mincarriers10_formeta.tsv.gz")
    sex <- 'maleonly'
  }else{
    mgb_file <- NA
  }  
  message("MGB file: ", mgb_file)
  
  # All three meta
  message("Running all three meta")
  cohort_files <- c(ukb_file, aou_file, mgb_file)
  cohort_files <- cohort_files[!is.na(cohort_files)]
  cohort_files_collapse <- paste0(cohort_files, collapse=" ")
  message(cohort_files_collapse)
  system(paste0("./UKBB_200KWES_CVD/metal_meta.sh UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_version ", cohort_files_collapse), intern=TRUE)
  if(length(cohort_files)>1){
    system(paste0("Rscript ./UKBB_200KWES_CVD/genotype_count_based_meta_analysis_dominant_underflow.R UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_version1.tbl ", 
                  "../UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv ",
                  0.1, " ", 0.05, " ", cohort_files_collapse), intern=TRUE) 
  }else{
    system(paste0("mv UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_version1.tbl ../UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv"), intern=TRUE)
    system(paste0("gzip ../UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv"))
  }
  
  # UKB AoU
  if(!(is.na(ukb_file) & is.na(aou_file))){
    message("Running UKB-AoU meta")
    if(is.na(mgb_file)){
      system(paste0("cp ../UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv.gz  ../UKB_AoUv8_meta_results_", code, "_", sex, "_adjusted.tsv.gz"), intern=TRUE)
    }else{
      cohort_files <- c(ukb_file, aou_file)
      cohort_files <- cohort_files[!is.na(cohort_files)]
      cohort_files_collapse <- paste0(cohort_files, collapse=" ")
      message(cohort_files_collapse)
      system(paste0("./UKBB_200KWES_CVD/metal_meta.sh UKB_AoUv8_meta_results_", code, "_", sex, "_version ", cohort_files_collapse), intern=TRUE)
      if(length(cohort_files)>1){
        system(paste0("Rscript ./UKBB_200KWES_CVD/genotype_count_based_meta_analysis_dominant_underflow.R UKB_AoUv8_meta_results_", code, "_", sex, "_version1.tbl ", 
                      "../UKB_AoUv8_meta_results_", code, "_", sex, "_adjusted.tsv ",
                      0.1, " ", 0.05, " ", cohort_files_collapse), intern=TRUE) 
      }else{
        system(paste0("mv UKB_AoUv8_meta_results_", code, "_", sex, "_version1.tbl ../UKB_AoUv8_meta_results_", code, "_", sex, "_adjusted.tsv"), intern=TRUE)
        system(paste0("gzip ../UKB_AoUv8_meta_results_", code, "_", sex, "_adjusted.tsv"))
      }
    }
  }
  
  # UKB MGB
  if(!(is.na(ukb_file) & is.na(mgb_file))){
    message("Running UKB-MGB meta")
    if(is.na(aou_file)){
      system(paste0("cp ../UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv.gz  ../UKB_MGB_meta_results_", code, "_", sex, "_adjusted.tsv.gz"), intern=TRUE)
    }else{
      cohort_files <- c(ukb_file, mgb_file)
      cohort_files <- cohort_files[!is.na(cohort_files)]
      cohort_files_collapse <- paste0(cohort_files, collapse=" ")
      message(cohort_files_collapse)
      system(paste0("./UKBB_200KWES_CVD/metal_meta.sh UKB_MGB_meta_results_", code, "_", sex, "_version ", cohort_files_collapse), intern=TRUE)
      if(length(cohort_files)>1){
        system(paste0("Rscript ./UKBB_200KWES_CVD/genotype_count_based_meta_analysis_dominant_underflow.R UKB_MGB_meta_results_", code, "_", sex, "_version1.tbl ", 
                      "../UKB_MGB_meta_results_", code, "_", sex, "_adjusted.tsv ",
                      0.1, " ", 0.05, " ", cohort_files_collapse), intern=TRUE) 
      }else{
        system(paste0("mv UKB_MGB_meta_results_", code, "_", sex, "_version1.tbl ../UKB_MGB_meta_results_", code, "_", sex, "_adjusted.tsv"), intern=TRUE)
        system(paste0("gzip ../UKB_MGB_meta_results_", code, "_", sex, "_adjusted.tsv"))
      }
    }
  }
  
  # AoU MGB
  if(!rerun){
    if(!(is.na(mgb_file) & is.na(aou_file))){
      message("Running AoU-MGB meta")
      if(is.na(ukb_file)){
        system(paste0("cp ../UKB_AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv.gz  ../AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv.gz"), intern=TRUE)
      }else{
        cohort_files <- c(aou_file, mgb_file)
        cohort_files <- cohort_files[!is.na(cohort_files)]
        cohort_files_collapse <- paste0(cohort_files, collapse=" ")
        message(cohort_files_collapse)
        system(paste0("./UKBB_200KWES_CVD/metal_meta.sh AoUv8_MGB_meta_results_", code, "_", sex, "_version ", cohort_files_collapse), intern=TRUE)
        if(length(cohort_files)>1){
          system(paste0("Rscript ./UKBB_200KWES_CVD/genotype_count_based_meta_analysis_dominant_underflow.R AoUv8_MGB_meta_results_", code, "_", sex, "_version1.tbl ", 
                        "../AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv ",
                        0.1, " ", 0.05, " ", cohort_files_collapse), intern=TRUE) 
        }else{
          system(paste0("mv AoUv8_MGB_meta_results_", code, "_", sex, "_version1.tbl ../AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv"), intern=TRUE)
          system(paste0("gzip ../AoUv8_MGB_meta_results_", code, "_", sex, "_adjusted.tsv"))
        }
      }
    }
  }
}
