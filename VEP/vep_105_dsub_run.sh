#!/bin/bash
set -o errexit

#python3 -m venv dsub_libs
#source dsub_libs/bin/activate
#pip install dsub

git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git
cd $HOME/TOPMed_AFib_pipeline/
git pull https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

# Author: Seung Hoan Choi <schoi@broadinstitute.org>
# Feb 7 2021
cd $HOME/TOPMed_AFib_pipeline/TOPMed_AFib_pipeline/VEP
gcloud config set project broad-ml4cvd

export PROJECT="$(gcloud config get-value project)"
#export GOOGLE_APPLICATION_CREDENTIALS="/Users/schoi/ellinor-ruff-bwh-timi-genomics-cb366b2d4166.json"

dsub \
   --project $PROJECT \
   --provider "google-cls-v2" \
   --use-private-address \
   --regions us-central1 us-east1 us-west1 \
   --boot-disk-size 500 \
   --disk-type pd-standard \
   --disk-size 1000 \
   --machine-type "n2-standard-8" \
   --image "gcr.io/broad-ml4cvd/vep:105_v6" \
   --skip \
   --logging "gs://ml4cvd/schoi/annotation/MGB_53K/log/" \
   --input VCF_FILE="gs://ml4cvd/schoi/annotation/MGB_53K/IBM_PHB_WES_callset_53K_Jan2022.filtered.0.vcf.gz" \
   --input LOFTEE_FILE="gs://ml4cvd/schoi/annotation/vep_105/GRCh38.tar" \
   --input VEP_FILE="gs://ml4cvd/schoi/annotation/vep_105/homo_sapiens_vep_105_GRCh38.tar.gz" \
   --input FASTA_FILE="gs://ml4cvd/schoi/annotation/vep_105/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.vep105.tar.gz" \
   --input PLUGIN_FILE="gs://ml4cvd/schoi/annotation/vep_105/loftee_Plugins.tar.gz" \
   --output OUTPUT_FILE="gs://ml4cvd/schoi/annotation/MGB_53K/annotated/IBM_PHB_WES_callset_53K_Jan2022.filtered.0.vcf.annotated.gz" \
   --output SEMAPHORE="gs://ml4cvd/schoi/annotation/MGB_53K/annotated/semaphore" \
   --script vep_105_dsub_script.sh \
   --wait


   To check the status, run:
     dstat --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220214-221911-31' --users 'schoi' --status '*'
   To cancel the job, run:
     ddel --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220214-221911-31' --users 'schoi'

     dstat --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-193044-79' --users 'schoi' --status '*'

     ddel --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-193044-79' --users 'schoi'

  dstat --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-163709-21' --users 'schoi' --status '*'

  ddel --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-163709-21' --users 'schoi'



  To check the status, run:
    dstat --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-164324-46' --users 'schoi' --status '*'
  To cancel the job, run:
    ddel --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-164324-46' --users 'schoi'

    To check the status, run:
  dstat --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-165028-21' --users 'schoi' --status '*'
To cancel the job, run:
  ddel --provider google-cls-v2 --project broad-ml4cvd --location us-central1 --jobs 'vep-105-ds--schoi--220209-165028-21' --users 'schoi'
