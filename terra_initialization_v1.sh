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
echo "BiocManager::install(c('SeqArray'),lib.loc ='/home/jupyter-user/.rpackages',  lib='/home/jupyter-user/.rpackages', dependencies=TRUE, clean=TRUE, INSTALL_opts='--no-docs --no-demo --byte-compile',version="3.8");" >> install.R
R CMD BATCH install.R
