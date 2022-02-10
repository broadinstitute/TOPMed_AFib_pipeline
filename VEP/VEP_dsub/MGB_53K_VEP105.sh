#gcloud compute ssh --zone=us-central1-a annot-shc --command="git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git;cd TOPMed_AFib_pipeline/VEP/VEP_dsub/;chmod 777 MGB_53K_VEP105.sh;./MGB_53K_VEP105.sh &"

docker pull gcr.io/google.com/cloudsdktool/cloud-sdk:latest

#### step 2
mkdir $HOME/vep_data
mkdir $HOME/vep_data/Plugins/
mkdir $HOME/vep_data/input/
mkdir $HOME/vep_data/output/


chmod a+rwx $HOME/vep_data
chmod a+rwx $HOME/vep_data/Plugins/
chmod a+rwx $HOME/vep_data/input/
chmod a+rwx $HOME/vep_data/output/
cd $HOME/vep_data/

#git clone -b grch38 https://github.com/konradjk/loftee.git
#cp -rf loftee/* Plugins/

#export PROJECT="$(gcloud config get-value project)"
#gsutil -u ellinor-ruff-bwh-timi-genomics -m cp gs://hail-us-vep/loftee-beta/GRCh38.tar gs://timi-cram/joint-callset/annotation/
#gsutil -u $PROJECT cat gs://hail-us-vep/loftee-beta/GRCh38.tar | tar -xf - -C $HOME/vep_data/
docker run -i -t -v $HOME/vep_data:/opt/vep/.vep gcr.io/google.com/cloudsdktool/cloud-sdk bash
export PROJECT="$(gcloud config get-value project)"

LOFTEE_FILE=gs://ml4cvd/schoi/annotation/vep_105/GRCh38.tar
VEP_FILE=gs://ml4cvd/schoi/annotation/vep_105/homo_sapiens_vep_105_GRCh38.tar.gz
FASTA_FILE=gs://ml4cvd/schoi/annotation/vep_105/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.vep105.tar.gz
PLUGIN_FILE=gs://ml4cvd/schoi/annotation/vep_105/loftee_Plugins.tar.gz
gsutil -m cp gs://ml4cvd/schoi/annotation/MGB_53K/IBM_PHB_WES_callset_53K_Jan2022.filtered.*.vcf.gz /opt/vep/.vep/
gsutil -m cat ${LOFTEE_FILE} | tar -xf - -C /opt/vep/.vep/
gsutil -m cat ${VEP_FILE} | tar -zxf - -C /opt/vep/.vep/
gsutil -m cat ${FASTA_FILE} | tar -zxf - -C /opt/vep/.vep/
gsutil -m cat ${PLUGIN_FILE} | tar -zxf - -C /opt/vep/.vep/Plugins/
exit

##### VEP annotation
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep gcr.io/broad-ml4cvd/vep:105_v5 bash
for num in {0..23}
do
/opt/vep/src/ensembl-vep/vep -i /opt/vep/.vep/input/IBM_PHB_WES_callset_53K_Jan2022.filtered.${num}.vcf.gz \
-o /opt/vep/.vep/output/IBM_PHB_WES_callset_53K_Jan2022.filtered.${num}.annotated.vcf.gz \
--format vcf \
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
--buffer_size 5000 \
--dir_plugins /opt/vep/.vep/Plugins \
--force_overwrite \
--fasta /opt/vep/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--plugin LoF,loftee_path:/opt/vep/.vep/Plugins,gerp_bigwig:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee.sql
done
exit

docker run -i -t -v $HOME/vep_data:/opt/vep/.vep gcr.io/google.com/cloudsdktool/cloud-sdk bash
export PROJECT="$(gcloud config get-value project)"
gsutil -m cp  $HOME/vep_data/output/IBM_PHB_WES_callset_53K_Jan2022.filtered.${num}.annotated.vcf.gz  gs://ml4cvd/schoi/annotation/MGB_53K/annotated/
exit
