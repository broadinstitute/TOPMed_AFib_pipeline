FROM ensemblorg/ensembl-vep:latest

MAINTAINER Seung Hoan Choi (schoi@broadinstitute.org)

RUN cpanm --local-lib=/opt/vep/perl5 DBD::SQLite::st::VirtualTable::PerlData
ENV PERL5LIB=$PERL5LIB://opt/vep/perl5/lib/perl5:${PERL5LIB}

USER root

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
    wget \
    make \
    gcc \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev

#
WORKDIR /usr/src

#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

RUN  wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
  tar jxf htslib-1.9.tar.bz2 && \
  rm htslib-1.9.tar.bz2 && \
  cd htslib-1.9 &&  \
  ./configure --prefix $(pwd) && \
	make

USER vep
ENV PATH=${PATH}:/usr/src/samtools-1.9
ENV PATH=${PATH}:/usr/src/htslib-1.9

WORKDIR /opt/vep/src/ensembl-vep
