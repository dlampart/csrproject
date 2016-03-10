#!/bin/bash

##################################
## Code/copyPreprocessedFiles.sh
##
## PreprocessedData/* /interimData/*
##################################
cp PreprocessedData/completeMotifFamilyTable.RData interimData/completeMotifFamilyTable.RData
cp PreprocessedData/motifFamilyTable.RData interimData/motifFamilyTable.RData
cp PreprocessedData/overallAveragedExpressionDirect.RDat interimData/overallAveragedExpressionDirect.RDat
cp PreprocessedData/unnormedMotifActivityDirect.RDat interimData/unnormedMotifActivityDirect.RDat

