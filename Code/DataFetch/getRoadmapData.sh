#!/bin/bash

cd Data
wget http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/57epigenomes.N.pc.gz
gunzip 57epigenomes.N.pc.gz
### remove universal reference sample E000
ar=(`head -1 57epigenomes.N.pc | perl -nle 's/\t/\n/mg && print' | perl -ne 'if($.>2){print}'`)
mkdir roadmapDHS
cd roadmapDHS
for id in ${ar[@]}
do
    echo $id
    wget http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/${id}-DNase.imputed.narrowPeak.bed.nPk.gz
done 
gunzip *
for id in ${ar[@]}
do
    echo $id
    sort-bed ${id}-DNase.imputed.narrowPeak.bed.nPk >  ${id}-DNase.imputed.narrowPeak.bed
    rm ${id}-DNase.imputed.narrowPeak.bed.nPk
done 
cd ../../

