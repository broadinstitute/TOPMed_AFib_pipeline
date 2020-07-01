########
######## KING relatedness
######## input enverionment
inputfile<-plinkfile

unrelated_prefix<-TOPMed_afib
unrelated_logfile<-outfile

duplicated_prefix<-TOPMed_afib
duplicated_logfile<-outfile

######
###### identify related individuals 2 degree
king -b ${inputfile}.bed --unrelated --degree 2 --prefix ${unrelated_prefix} > ${unrelated_logfile}

######
###### duplicated individuals
king -b ${inputfile}.bed --duplicate --prefix ${duplicated_prefix} > ${duplicated_logfile}

########
########
######## plink kinship relatedness
######## input enverionment
kinship_prefix<-plink_kinship
unrelated_logfile<-outfile

######
###### kinship matrix
plink --bfile ${inputfile} --make-rel square --out {kinship_prefix} > ${unrelated_logfile}

########
########
######## PC calculation
######## input enverionment
unrelated_output<-plinkname
related_output<-plinkname
tg_input<-plinkname
unrelated_tg_output<-plinkname

pc_step1_prefix<-prefix
pc_step1_log<-logfile


###### create unrelated data set
plink --bfile ${inputfile} --keep ${unrelated_prefix}unrelated.txt --keep-allele-order --make-bed --out ${unrelated_output}

#### create related data set
plink --bfile ${inputfile} --remove ${unrelated_prefix}unrelated.txt --keep-allele-order --make-bed --out ${related_output}

##### combine unrelated TOPMed data set and 1000G data set
plink --bfile ${unrelated_output} --bmerge ${tg_input} --keep-allele-order --make-bed --out ${unrelated_tg_output}


##### PC step 1
##### PC analysis using unrelated ( + 1000G) individual

flashpca --bfile ${unrelated_tg_output} --ndim 20 \
--outload ${pc_step1_prefix}_unrelated_loadings.txt \
--outmeansd ${pc_step1_prefix}_unrelated_meansd.txt \
--outpc ${pc_step1_prefix}_unrelated_pcs.txt \
--outpve ${pc_step1_prefix}_unrelated_pve.txt \
--outvec ${pc_step1_prefix}_unrelated_eigenvalues.txt \
--outval ${pc_step1_prefix}_unrelated_eigenvectors.txt > ${pc_step1_log}


##### PC step 2
##### PC analysis projection to related individuals

flashpca --bfile ${related_output} --ndim 20 \
--project --inmeansd ${pc_step1_prefix}_unrelated_meansd.txt \
--inload ${pc_step1_prefix}_unrelated_loadings.txt \
--outproj ${pc_step1_prefix}_unrelated_projections.txt -v > ${pc_step1_prefix}_step2_projection.out


####
#### version 2 TOPMed only

pc_step1_prefix2<-prefix
pc_step1_log2<-outfile


flashpca --bfile ${unrelated_output} --ndim 20 \
--outload ${pc_step1_prefix2}_loadings.txt \
--outmeansd ${pc_step1_prefix2}_meansd.txt \
--outpc ${pc_step1_prefix2}_pcs.txt \
--outpve ${pc_step1_prefix2}_pve.txt \
--outvec ${pc_step1_prefix2}_eigenvalues.txt \
--outval ${pc_step1_prefix2}_eigenvectors.txt > ${pc_step1_log2}


flashpca --bfile ${related_output} --ndim 20 \
--project --inmeansd ${pc_step1_prefix2}_unrelated_meansd.txt \
--inload ${pc_step1_prefix2}_unrelated_loadings.txt \
--outproj ${pc_step1_prefix2}_unrelated_projections.txt -v > ${pc_step1_prefix2}_step2_projection.out



echo "=========================================================="
echo "Finished on       : $(date)"
echo "=========================================================="
