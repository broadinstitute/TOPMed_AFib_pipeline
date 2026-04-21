chunk_num=$1
n_chunks=$2

mkdir tmp/
cd tmp/
  
git clone https://github.com/statgen/METAL
cd METAL
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install
cd ../..

git clone https://github.com/seanjosephjurgens/UKBB_200KWES_CVD/
chmod +x ./UKBB_200KWES_CVD/metal_meta.sh
chmod +x ./UKBB_200KWES_CVD/genotype_count_based_meta_analysis_dominant.R

git clone --branch patch-1 https://github.com/seanjosephjurgens/TOPMed_AFib_pipeline/
chmod +x ./TOPMed_AFib_pipeline/DNANexus/PheWAS/GC_based_meta/GC_based_meta_UKB_AoUv8_MGB.R

dx download exome-seq:/sjj/projects/phewas/v2/sum_stats_collected/UKB/for_meta/*formeta.tar.gz
tar -xzf UKB_sumstats_formeta.tar.gz
dx download exome-seq:/sjj/projects/phewas/v2/sum_stats_collected/AoUv8/for_meta/*formeta.tar.gz
tar -xzf AoUv8_sumstats_formeta.tar.gz
dx download exome-seq:/sjj/projects/phewas/v2/sum_stats_collected/MGB/for_meta/*formeta.tar.gz
tar -xzf MGB_sumstats_formeta.tar.gz

Rscript ./TOPMed_AFib_pipeline/DNANexus/PheWAS/GC_based_meta/GC_based_meta_UKB_AoUv8_MGB.R ${chunk_num} ${n_chunks}

cd ..
rm -rf tmp/
rm *.sh

#dx upload * --destination exome-seq:/sjj/projects/phewas/v2/sum_stats_collected/UKB_AoUv8_MGB_META/
