#! /bin/bash
set -o errexit

# Author: Seung Hoan Choi <schoi@broadinstitute.org>
# Feb 07 2022

mkdir $HOME/vep_data
mkdir $HOME/vep_data/Plugins/

chmod a+rwx $HOME/vep_data
chmod a+rwx $HOME/vep_data/Plugins/

tar -xf ${MYBUCKET}/schoi/annotation/vep_105/GRCh38.tar -C /opt/vep/.vep/
tar -zxf ${MYBUCKET}/schoi/annotation/vep_105/homo_sapiens_vep_105_GRCh38.tar.gz -C /opt/vep/.vep/
tar -zxf ${MYBUCKET}/schoi/annotation/vep_105/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.vep105.tar.gz -C /opt/vep/.vep/
tar -zxf ${MYBUCKET}/schoi/annotation/vep_105/loftee_Plugins.tar.gz -C /opt/vep/.vep/Plugins/

#### instlal packages
cpanm DBD::SQLite::VirtualTable::PerlData
export PERL5LIB=$PERL5LIB://opt/vep/perl5/lib/perl5:${PERL5LIB}


# TODO: Consider setting --lowmem if #traits > 10 as per
# https://rgcgithub.github.io/regenie/options/#input

# For bsize choice, please see https://rgcgithub.github.io/regenie/faq/ for
# rationale.
/opt/vep/src/ensembl-vep/vep \
-i ${VCF_FILE} \
-o ${OUT_FILE} \
--compress_output gzip \
--assembly GRCh38 --species homo_sapiens \
--offline --cache \
--no_stats \
--minimal \
--everything \
--allele_number \
--show_ref_allele \
--af_gnomad \
--canonical --tab \
--dir_plugins /opt/vep/.vep/Plugins \
--force_overwrite \
--fasta /opt/vep/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--plugin LoF,loftee_path:/opt/vep/.vep/Plugins,gerp_bigwig:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee.sql


echo "OK" > ${SEMAPHORE}
