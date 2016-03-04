#!/bin/bash
cd Data/CELs
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE15nnn/GSE15805/suppl/GSE15805_RAW.tar
tar xvf GSE15805_RAW.tar 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19090/suppl/GSE19090_RAW.tar
tar xvf GSE19090_RAW.tar
gunzip *.gz
cd ..
wget http://www.stemformatics.org/mappings/mapping_71.txt
