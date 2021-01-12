from quay.io/biocontainers/bioconductor-seqarray:1.30.0--r40h5f743cb_0

MAINTAINER Seung Hoan Choi (schoi@broadinstitute.org)

RUN apt-get update
RUN apt-get -y install git

RUN git clone https://github.com/seuchoi/vcf2gds.git && cd ./vcf2gds && git pull origin master
