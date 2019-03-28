# Base Image
FROM biocontainers/bcftools:v1.5_cv2

# Maintainer
MAINTAINER Seung Hoan Choi <seuchoi@gmail.com>

RUN git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git && cd ./TOPMed_AFib_pipeline && git pull origin master
