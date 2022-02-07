#!/bin/bash
set -o errexit

#python3 -m venv dsub_libs
#source dsub_libs/bin/activate

# Author: Seung Hoan Choi <schoi@broadinstitute.org>
# Feb 7 2021
cd /Users/schoi/github/TOPMed_AFib_pipeline/VEP/
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
   --image "gcr.io/ellinor-ruff-bwh-timi-genomics/vep:105" \
   --skip \
   --wait \
   --mount MYBUCKET=gs://ml4cvd \
   --logging "gs://ml4cvd/schoi/annotation/vep_105/" \
   --input VCF_FILE="gs://timi-cram/joint-callset/annotation/MGB_53K/IBM_PHB_WES_callset_53K_Jan2022.filtered.0.vcf.gz" \
   --output OUTPUT_FILE="gs://ml4cvd/schoi/annotation/vep_105/annotated/IBM_PHB_WES_callset_53K_Jan2022.filtered.0.vcf.annotated.gz" \
   --output SEMAPHORE="gs://ml4cvd/schoi/annotation/vep_105/annotated/semaphore" \
   --script vep_105_dsub_script.sh
