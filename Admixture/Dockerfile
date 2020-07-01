# Base Image
FROM rtibiocloud/plink:v1.9_bc8fbd1

# Maintainer
MAINTAINER Seung Hoan Choi <seuchoi@gmail.com>

## pulling my file
RUN git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git && cd ./TOPMed_AFib_pipeline && git pull origin master

##
RUN ADMIXTURE_VERSION=1.3.0 && \
    ADMIXTURE_GZ_PATH=/tmp/admixture-$ADMIXTURE_VERSION.tar.gz && \
    curl -L -o $ADMIXTURE_GZ_PATH http://software.genetics.ucla.edu/admixture/binaries/admixture_linux-1.3.0.tar.gz && \
    mkdir /tmp/admixture/ && \
    tar -xvzf $ADMIXTURE_GZ_PATH -C /tmp/admixture/ && \
    cp /tmp/admixture/admixture_linux-1.3.0/* ../ && \
    ln -s /tmp/admixture/admixture /bin/admixture 
