#!/bin/bash

#mkdir Data/genomicSequence
cd Data/genomicSequence
for j in `seq 1 1 22`
do
    echo $j    
    wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${j}.fa.gz
    gunzip chr${j}.fa.gz
done