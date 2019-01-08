#!/bin/sh

#update java
sudo apt install openjdk-8-jdk

#install fastq-join
git clone https://github.com/brwnj/fastq-join
cd fastq-join
make
sudo mv fastq-join /usr/local/bin/
cd ..

#install bwa 
sudo apt-get install bwa 




#install CLARK
wget http://clark.cs.ucr.edu/Download/CLARKV1.2.5.tar.gz
tar -xzvf CLARKV1.2.5.tar.gz 
cd CLARKSCV1.2.5.1/
./install.sh 
./set_targets.sh  DIR_DB bacteria viruses fungi human
./updateTaxonomy.sh 
gunzip ../example/sample1_meta_R* 
./classify_metagenome.sh --light -O ../example/setup/sample1_meta_R1.fastq -R ../example/setup/sample1_meta_R1
rm -r DIR_DB/Bacteria/
rm -r DIR_DB/Viruses/
rm -r DIR_DB/Fungi/
rm -r DIR_DB/Human/
cd taxonomy 
rm !(nodes.dmp)
cd ../..


#install krona tools
wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar -vxf KronaTools-2.7.tar
cd KronaTools-2.7/
sudo ./install.pl
cd .. 

#download human genome 
mkdir host
cd host
wget         ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
cat *.gz > human.fasta.gz
rm !(human.fasta.gz)
bwa index human.gz
cd ..


#install nextflow
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# install fastqc 
brew update
brew install fastqc

pip install --upgrade pip
hash -d pip
pip install --user multiqc

#install seqtk
sudo apt install seqtk

#?symlink scripts/clark etc 

