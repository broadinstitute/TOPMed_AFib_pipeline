#!/usr/bin/env Rscript

#### binary traits
args=(commandArgs(TRUE))
study_result_file_string=as.character(args[1])
study_result_name_string=as.character(args[2])
min_carrier_study=as.numeric(args[3])
min_carrier_meta=as.numeric(args[4])
chr=as.numeric(args[5])

cat('\nReading in packages for analysis...\n')
### library(GENESIS)
### library(data.table)
### library(SeqArray)
### library(SeqVarTools)

.libPaths(c("rpackages4_1_3",.libPaths()))

#git clone --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#git pull --branch v1.2 https://github.com/seanjosephjurgens/UKBB_200KWES_CVD.git
#source("/medpop/afib/sjurgens/Rscripts/association_source_v2.R")

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("UKBB_200KWES_CVD/score_based_meta_analysis_methods.R")
source("UKBB_200KWES_CVD/Cauchy_test.R")

study_file_split <- strsplit(study_result_file_string, ";")[[1]]
study_name_split <- strsplit(study_result_name_string, ";")[[1]]

study_list <- list()
cat('Working on results for chr', chr, '...\n')
for(i in c(1:length(study_file_split))){
    inter <- fread(study_file_split[i], stringsAsFactors=F, data.table=F)
    inter <- inter[inter$chr == chr, ]
    total <- paste0(inter$Phenotype, "__", inter$gene, "__", inter$mask)
    rm <- which(duplicated(total))
    if(length(rm)>0){inter <- inter[-rm, ]; total <- total[-rm]}
    rownames(inter) <- total
    cat('\t', paste0(study_name_split[i]), '\n')
    head(inter)
    study_list[[paste0(study_name_split[i])]] <- inter
}

score_col_vec <- c(rep('Score', length(study_file_split)))
pval_col_vec <- c(rep('SPA.pval', length(study_file_split)))
score.se_col_vec <- c(rep('Score.SE', length(study_file_split)))
est_col_vec <- c(rep('Est', length(study_file_split)))
est.se_col_vec <- c(rep('Est.SE', length(study_file_split)))
est_type="logistic"
recalc_variance_vec <- c(rep(T, length(study_file_split)))
mincarriers_col_vec <- c(rep('cMAC', length(study_file_split)))
mincarriers_num_vec <- c(rep(min_carrier_study, length(study_file_split)))
meta_mincarriers_num=min_carrier_meta

study_result_name_string <- gsub(";", "_", study_result_name_string)
outfile <- paste0(study_result_name_string, "__10CarriersPerStudy_20Overall_all_tests_chr", chr, '.tsv')
cat('Writing results to ', outfile, '...\n')
cat('Starting meta-analysis...\n')
meta <- score_meta(single_variant=F,
                   study_summary_data_list=study_list,
                   score_col_vec=score_col_vec,
                   pval_col_vec=pval_col_vec,
                   score.se_col_vec=score.se_col_vec,
                   est_col_vec=est_col_vec,
                   est.se_col_vec=est.se_col_vec,
                   est_type="logistic",
                   calc_raw_meta_fixedeffects_odds_ratio=T,
                   recalc_variance_vec=recalc_variance_vec,
                   mincarriers_col_vec=mincarriers_col_vec,
                   mincarriers_num_vec=mincarriers_num_vec,
                   meta_mincarriers_num=meta_mincarriers_num,
                   save_output=F,
                   outfile=outfile,
                   make_figures=F,
                   min_number_of_studies_contributing=1,
                   n.cores=1
)

sessionInfo()
quit("no")
