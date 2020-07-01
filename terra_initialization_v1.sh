#!/bin/bash
set -e
### setup envirnment
### setup htslib
cd /tmp/
git clone git://github.com/samtools/htslib.git
cd htslib
make
ln -s /tmp/htslib/bgzip /bin/bgzip
ln -s /tmp/htslib/htsfile /bin/htsfile
ln -s /tmp/htslib/tabix /bin/tabix

### setup bcftools
cd /tmp/
git clone git://github.com/samtools/bcftools.git
cd bcftools
make
ln -s /tmp/bcftools/bcftools /bin/bcftools

#### setup plink
PLINK_VERSION=1.9.0Beta
PLINK_ZIP_PATH=/tmp/plink-$PLINK_VERSION.zip
curl -L -o $PLINK_ZIP_PATH http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip
unzip -o $PLINK_ZIP_PATH -d /tmp/plink/
ln -s /tmp/plink/plink /bin/plink


#### setup library
#### install SeqArray
cd /tmp/
echo 'install.packages("BiocManager",lib="/home/jupyter-user/.rpackages")' > install.R
#echo "BiocManager::install(c('SeqArray'),lib.loc ='/home/jupyter-user/.rpackages',  lib='/home/jupyter-user/.rpackages', dependencies=TRUE, clean=TRUE, INSTALL_opts='--no-docs --no-demo --byte-compile',version="3.10");" >> install.R
echo "BiocManager::install(c('SeqArray'),lib.loc ='/home/jupyter-user/.rpackages',  lib='/home/jupyter-user/.rpackages', dependencies=TRUE, clean=TRUE, INSTALL_opts='--no-docs --no-demo --byte-compile');" >> install.R
R CMD BATCH install.R

##### set up the pipeline
cd /home/jupyter-user/
git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git


#####
##### plink 2 install
PLINK_VERSION=2.0.Alpha
PLINK_ZIP_PATH=/tmp/plink-$PLINK_VERSION.zip
curl -L -o $PLINK_ZIP_PATH http://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20191104.zip
mkdir /tmp/plink2/
unzip -o $PLINK_ZIP_PATH -d /tmp/plink2/
ln -s /tmp/plink2/plink2 /bin/plink2

######
###### admixture install
ADMIXTURE_VERSION=1.3.0
ADMIXTURE_GZ_PATH=/tmp/admixture-$ADMIXTURE_VERSION.tar.gz

curl -L -o $ADMIXTURE_GZ_PATH http://software.genetics.ucla.edu/admixture/binaries/admixture_linux-1.3.0.tar.gz
mkdir /tmp/admixture/
tar -xvzf $ADMIXTURE_GZ_PATH -C /tmp/admixture/
cp /tmp/admixture/admixture_linux-1.3.0/* ../
ln -s /tmp/admixture/admixture /bin/admixture


######
###### install king
KING_VERSION=2.2
KING_GZ_PATH=/tmp/king-$KING_VERSION.tar.gz

curl -L -o $KING_GZ_PATH http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz
mkdir /tmp/king/
tar -xvzf $KING_GZ_PATH -C /tmp/king/
ln -s /tmp/admixture/admixture /bin/admixture
