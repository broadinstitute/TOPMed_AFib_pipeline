#!/bin/bash
set -o errexit

#python3 -m venv dsub_libs
#source dsub_libs/bin/activate
#pip install dsub

git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git
# Author: Seung Hoan Choi <schoi@broadinstitute.org>
# Feb 7 2021
cd TOPMed_AFib_pipeline/VEP
export PROJECT="$(gcloud config get-value project)"
#export GOOGLE_APPLICATION_CREDENTIALS="/Users/schoi/ellinor-ruff-bwh-timi-genomics-cb366b2d4166.json"

dsub \
   --project $PROJECT \
   --provider "google-cls-v2" \
   --use-private-address \
   --regions us-central1 us-east1 us-west1 \
   --disk-type pd-ssd \
   --disk-size 1000 \
   --machine-type "n2-custom-80-40960" \
   --image "gcr.io/broad-ml4cvd/vep:105" \
   --skip \
   --wait \
   --logging "gs://ml4cvd/schoi/annotation/MGB_53K/log/" \
   --input VCF_FILE="gs://ml4cvd/schoi/annotation/MGB_53K/IBM_PHB_WES_callset_53K_Jan2022.filtered.0.vcf.gz" \
   --input LOFTEE_FILE="gs://ml4cvd/schoi/annotation/vep_105/GRCh38.tar" \
   --input VEP_FILE="gs://ml4cvd/schoi/annotation/vep_105/homo_sapiens_vep_105_GRCh38.tar.gz" \
   --input FASTA_FILE="gs://ml4cvd/schoi/annotation/vep_105/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.vep105.tar.gz" \
   --input PLUGIN_FILE="gs://ml4cvd/schoi/annotation/vep_105/loftee_Plugins.tar.gz" \
   --output OUTPUT_FILE="gs://ml4cvd/schoi/annotation/MGB_53K/annotated/IBM_PHB_WES_callset_53K_Jan2022.filtered.0.vcf.annotated.gz" \
   --output SEMAPHORE="gs://ml4cvd/schoi/annotation/MGB_53K/annotated/semaphore" \
   --script vep_105_dsub_script.sh
