#!/bin/bash

#install homebrew 

URL_BREW='https://raw.githubusercontent.com/Homebrew/install/master/install'

echo -n '- Installing brew ... '
echo | /usr/bin/ruby -e "$(curl -fsSL $URL_BREW)" > /dev/null
if [ $? -eq 0 ]; then echo 'OK'; else echo 'NG'; fi

#update java
brew cask install java



#install bwa 
brew install bwa


sudo easy_install pip

pip install futures

#install CLARK
pip install  --user futures

curl -O  http://clark.cs.ucr.edu/Download/CLARKV1.2.5.tar.gz
tar -xzvf CLARKV1.2.5.tar.gz 
cd CLARKSCV*
./install.sh
mkdir DIR_DB
cd DIR_DB 
cp ../../resources/xa* .
cat xa* > C-VRDB.tar.gz
tar -xzvf C-VRDB.tar.gz

cd ../

./set_targets.sh  DIR_DB bacteria fungi human custom
./updateTaxonomy.sh 
gunzip ../example/sample1_meta_R* 
./classify_metagenome.sh --light -O ../example/setup/1_S1_L001_R1_001.fastq -R ../example/setup/sample1_meta_R1
rm -r DIR_DB/Bacteria/
rm -r DIR_DB/Custom/
rm -r DIR_DB/Fungi/
rm -r DIR_DB/Human/
cd DIR_DB/taxonomy 
shopt -s extglob 
rm -- !(nodes.dmp)
cd ../../..


#install krona tools
curl -O  https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar -vxf KronaTools-2.7.tar
cd KronaTools-2.7/
sudo ./install.pl
./updateTaxonomy.sh 
cd .. 

#download human genome 
mkdir host
cd host
curl -O  ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
cat *.gz > human.fasta.gz
rm !(human.fasta.gz)
bwa index human.fasta.gz
cd ..


#install nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# install fastqc 
brew update
brew install fastqc

pip install --upgrade pip
hash -d pip
pip install --user multiqc

#install seqtk
brew install seqtk

echo "Installation finish $(date)" 

