#!/bin/bash

cd Data
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.N.pc.gz
gunzip 57epigenomes.N.pc.gz
### remove universal reference sample E000
ar=(`head -1 57epigenomes.N.pc | perl -nle 's/\t/\n/mg && print' | perl -ne 'if($.>2){print}'`)
mkdir imputedDHS
cd imputedDHS
for id in ${ar[@]}
do
    echo $id
    wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/${id}-DNase.imputed.narrowPeak.bed.gz
#wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/${id}-DNase.imputed.narrowPeak.bed
done 
gunzip *
