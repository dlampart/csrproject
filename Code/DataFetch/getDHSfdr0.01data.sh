#!/bin/bash

#######################################################
## Code/DataFetch/getDHSfdr0.01data.sh
## input:
## output:Data/DHSfdr0.01/wgEncode(Uw|Duke)Dnase(\w+).fdr01peaks.hg19.bb.bed
## output:Data/multi-tissue.master.ntypes.simple.hg19.bed

## explanation: Downloads DHS data
## additionally, download Data/multi-tissue.master.ntypes.simple.hg19.bed
#######################################################
cd Data
wget ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/combined_peaks/multi-tissue.master.ntypes.simple.hg19.bed
more multi-tissue.master.ntypes.simple.hg19.bed | perl -pe 's/MCV\-/.\t/' > ../interimData/nrOfCelllines.bed
cd ..
#mkdir Data/DHSfdr0.01
cd Data/DHSfdr0.01
touch index.html
rm index.html
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/fdrPeaks/
my_names=(`more index.html | perl -ne '/(wgEncode.+?bb)/ && print "$1\n"'`)
chmod +x ../../ExternalCode/bigBedToBed
for my_name in ${my_names[@]}
do
echo $my_name
wget http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/fdrPeaks/$my_name
../../ExternalCode/bigBedToBed $my_name $my_name.bed
done
rm *.bb
