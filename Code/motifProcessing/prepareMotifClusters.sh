#!/bin/bash                                                                                      
#######################################################
### Code/motifProcessing/prepareMotifClusters.sh
###
### input: http://tfclass.bioinf.med.uni-goettingen.de/suplementary/TFClass.obo
### output interimData/tfClassifications.txt
### columns owl-id; uniprot-id
### 
### Details: 
### TFClass.obo contains for each transcription factor its tf-subfamily and family membership.
### Since tfs are organized within a tree (almost), each node has  one parent and the whole 
### membership structure is encoded in the tf-id. 
### so we only need to keep the lowest level tf-ids 
#######################################################
## extract all subfamily levels and ENSG                                       
cd Data
wget http://tfclass.bioinf.med.uni-goettingen.de/suplementary/TFClass.obo
cd ../
### keep only ids that represent genes and the corresponding uniprot id
more Data/TFClass.obo | perl -ne '/id: \d+\.\d+\.\d+\.\d+\.\d+|xref: UNIPROT/ && print' > interimData/tmp
### get both ids onto same line
more interimData/tmp | perl -pe 's/(\.\d+\.\d+)\n/$1\t/mg' > interimData/tmp2
### remove alternative ids
more interimData/tmp2 | perl -pe 's/alt_id.+\t//' > interimData/tmp3
### remove id signifier strings:
more interimData/tmp3 | perl -pe 's/id: //' | perl -pe 's/xref: UNIPROT://' >  interimData/tfClassifications.txt
### 
