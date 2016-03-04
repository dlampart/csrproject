#!/bin/bash

## input: Data/mapping_71.txt (Downloaded from http://www.stemformatics.org/contents/download_mappings)
## input: Data/martEnsg2hgnc.txt (downloaded from biomart)
## input: Data/UniProtGeneSymbolTable.tbl
##  (downloaded from http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22, columns: Entry; Entry name; Gene names(primary))


###################################
## this script calls the 
## the processing script used in
## the paper 'papername'. 
##
## this pipeline uses precomputed data 
## matrices in folder PreprocessedData/ .
## Calling  Code/pipelineWithoutPreprocessing.R 
## will overwrite some processed data in interimData/
## with precomputed data in PreprocessedData/ .
###################################
##
###################################
bash Code/preparePathStructures.sh


##################################
## Code/copyPreprocessedFiles.sh
##
## PreprocessedData/* /interimData/*
##################################
bash Code/copyPreprocessedFiles.sh

 



#######################################################
## Code/prepareNormalizedMatricesDirectFinal.R
##
## uses Code/normalizationFunctions.R
## input: interimData/unnormedMotifActivityDirect.RDat
## input: interimData/overallAveragedExpressionDirect.RDat
## output: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output: interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat
##
## explanations: normalizes both matrices and removes principal compnent with
## largest eigenvector from the motif matrix.
#######################################################
Rscript Code/prepareNormalizedMatricesDirectFinal.R

#######################################################
## Code/makeFigures/runAllMethodsOnlyForTFsFinal.R
##
## uses: Code/makeFigures/helperFunctionsFinal.R
## uses: Code/fast_lmm.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output: interimData/allResMinus.RDat
## output: interimData/allRes0.RDat
## output: interimData/allRes1.RDat
##
## explanations: runs fast_lmm only for Transcription factors
#######################################################
Rscript Code/makeFigures/runAllMethodsOnlyForTFsFinal.R

#######################################################
## Code/makeFigures/prepareRankingTablesFinal.R
##
## uses: Code/makeFigures/helperFunctionsFinal.R
## interimData/completeMotifFamilyTable.RData
## input: interimData/allResMinus.RDat
## input: interimData/allRes0.RDat
## input: interimData/allRes1.RDat
## output: interimData/alldfMinus.RDat
## output: interimData/alldf0.RDat
## output: interimData/alldf1.RDat
##
## explanation: prepares tables motif rank score
#######################################################
Rscript Code/makeFigures/prepareRankingTablesFinal.R

#######################################################
## Code/makeFigures/runAllScoresFinal.R
##
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output: interimData/completeRun.RDat
## explanation: runs fast_lmm across all motifs and expression values
#######################################################
Rscript Code/makeFigures/runAllScoresFinal.R

#######################################################
## Code/makeFigures/prepareRankingTablesPerSubFamilyFinal.R
##
## use: Code/makeFigures/helperFunctionsFinal.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/allRes1.RDat
## input: interimData/allRes0.RDat
## input: interimData/allResMinus.RDat
## input: interimData/alldf1.RDat
## input: interimData/alldf0.RDat
## input: interimData/alldfMinus.RDat
## output: interimData/alldf0SubFam.RDat
## output: interimData/alldf1SubFam.RDat
## output: interimData/alldfMinusSubFam.RDat
#######################################################
Rscript Code/makeFigures/prepareRankingTablesPerSubFamilyFinal.R

#######################################################
## Code/makeFigures/addVarianceEffectToSubFamilyTable.R
##
## use: Code/makeFigures/helperFunctionsFinal.R
## input: prepareRankingTablesPerSubFamilyFinal.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/alldf1SubFam.RDat
## output: interimData/alldf1SubFamWithVariance.RDat
#######################################################
Rscript Code/makeFigures/addVarianceEffectToSubFamilyTable.R

#######################################################
## Code/makeFigures/makeFigureCheckInflationControlFinal.R
##
## uses Code/naive_regression.R
## uses Code/fast_lmm.R
## input: interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output: PaperDocs/Images/figControlCheckPanel(A|B|C)_COE1.svg (Fig 2)
#######################################################
Rscript Code/makeFigures/makeFigureCheckInflationControlFinal.R

#######################################################
## Code/makeFigures/showPowerIncreasePerSubfamilyFinal.R
##
## input: interimData/alldfMinusSubFam.RDat
## input: interimData/alldf0SubFam.RDat
## input: interimData/alldf1SubFam.RDat
## output: PaperDocs/Images/showPowerCumulativePerSubfamily.svg (Fig 4)
#######################################################
Rscript Code/makeFigures/showPowerIncreasePerSubfamilyFinal.R

#######################################################
## Code/makeFigures/preparePioneerPlotFinal.R
## uses Code/makeFigures/helperFunctionsFinal.R
## interimData/completeMotifFamilyTable.RDatai
## interimData/alldf1SubFamWithVariance.RDat
## output: PaperDocs/Images/pioneerEnrichment.svg (Fig 5)
#######################################################
Rscript Code/makeFigures/preparePioneerPlotFinal.R

#######################################################
## Code/makeFigures/makeFigurePlotQQ_GRlikeFinal.R
##
## uses: Code/fast_lmm.R
## input:interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output:PaperDocs/Images/qqPlotER-Like_GCR.svg (Fig 6)
#######################################################
Rscript Code/makeFigures/makeFigurePlotQQ_GRlikeFinal.R

#######################################################
## Code/makeFigures/makeFigureDisplayPOU5F1.R
##
## uses: Code/fast_lmm.R
## input:interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output:PaperDocs/Images/figure_PO5F.svg (SupFig 2)
######################################################
Rscript Code/makeFigures/makeFigureDisplayPOU5F1.R

#######################################################
## Code/makeFigures/showCorrelationPlots.R
##
## input: interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat
## output:PaperDocs/Images/cor2Motif0.svg
## output:PaperDocs/Images/cor2Motif1.svg
## output:PaperDocs/Images/cor2Expr.svg
## output: PaperDocs/Images/eigenValues.svg
## (supFig 4)
######################################################
Rscript Code/makeFigures/showCorrelationPlots.R

#######################################################
## Code/makeFigures/makeVarianceEffectPlot.R
##
## input: interimData/alldf1.RDat
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/overallAveragedExpressionDirect.RDat
## input: interimData/alldf1SubFamWithVariance.RDat
## output: PaperDocs/Images/showVarianceEffect.pdf
## output: PaperDocs/Images/showVarianceEffectSubfamily.pdf
## (supFig 3)
######################################################
Rscript Code/makeFigures/makeVarianceEffectPlot.R
