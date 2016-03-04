#!/bin/bash
cd Data/chipSeq
cd Haib
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/
ar=(`more index.html | grep href | grep broadPeak | perl -ne '/(wg.+?gz)/ && print "$1\n"'`)
for el in "${ar[@]}"
do
    echo $el
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/${el}
done

cd ..
cd Sydh
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/
ar=(`more index.html | grep href | grep narrowPeak | perl -ne '/(wg.+?gz)/ && print "$1\n"'`)
for el in "${ar[@]}"
do
    echo $el
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/${el}
done
cd ..
wget https://raw.githubusercontent.com/ENCODE-DCC/dccMetadataImport/master/data/ucsc/mdb/human/wgEncodeHaibTfbs.ra
wget https://raw.githubusercontent.com/ENCODE-DCC/dccMetadataImport/master/data/ucsc/mdb/human/wgEncodeSydhTfbs.ra

more wgEncodeHaibTfbs.ra | perl -ne '!/^$/ && print' | perl -nle 'if(/^metaObject/){$c=$c+1}{print "$c\t$_";}' | perl -ple 's/ /\t/' | cut -f1,2,3 > wgEncodeHaibTfbsPrepared.txt

more wgEncodeSydhTfbs.ra | perl -ne '!/^$/ && print' | perl -nle 'if(/^metaObject/){$c=$c+1}{print "$c\t$_";}' | perl -ple 's/ /\t/' | cut -f1,2,3 > wgEncodeSydhTfbsPrepared.txt


