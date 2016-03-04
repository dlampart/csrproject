#!/bin/bash

if [ `ls | grep PaperDocs | wc -l` -eq "0" ];then
###
###
    echo "set up ExternalCode/ paths"
    mkdir ExternalCode/
###
###
    echo "set up Data/ paths"
    mkdir Data/motifData/
    mkdir Data/DHSfdr0.01/
    mkdir Data/genomicSequence/
    mkdir Data/CELs/
    mkdir Data/chipSeq/
    mkdir Data/chipSeq/Haib/
    mkdir Data/chipSeq/Sydh/
###
###
    echo "set up interimData/ paths"
    mkdir interimData/
    mkdir interimData/DHSnormalized/
    mkdir interimData/OverlappedCELs/
    mkdir interimData/motifInstances/
    mkdir interimData/motifInstances/HOCOMOCOSPLITTED/
    mkdir interimData/motifInstances/HOCOMOCOBEDS/
    mkdir interimData/motifInstances/HOCOMOCOHIGHCONF/
###
###
    echo "set up PaperDocs/ paths"
    mkdir PaperDocs/
    mkdir PaperDocs/Images/
fi