# Base Image
FROM avelior/plink2:latest

# Maintainer
MAINTAINER Seung Hoan Choi <seuchoi@gmail.com>

## pulling my file
RUN git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git && cd ./TOPMed_AFib_pipeline && git pull origin master

##
RUN wget http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz && \
    tar -xzvf Linux-king.tar.gz && \
    cp king /bin/king
