#!/bin/bash

###################################
## this script calls the 
## the processing script used in
## the paper
## 'Genome-wide association between transcription factor expression and chromatin accessibility reveals chromatin state regulators'
## in order.
##
## while it should allow to replicate 
## results with relative ease it's not guaranteed to work out-of-the-box.
## 
## notworthy failure points are installation
## of R packages, external software
## and download of external data.
##
## Be warned: some of the steps are
## computationally very intensive.
## A second pipeline (script Code/pipelineWithoutPreprocessing.R)
## is available that skips preprocessing. It uses
## precomputed data matrix in folder PreprocessedData/ .
## Calling  Code/pipelineWithoutPreprocessing.R 
## will overwrite processed data in interimData/
## with processed data in PreprocessedData/ .
###################################

## Data not automatically downloaded:
## input: Data/mapping_71.txt (Downloaded from http://www.stemformatics.org/contents/download_mappings)
## input: Data/martEnsg2hgnc.txt (downloaded from biomart)
## input: Data/UniProtGeneSymbolTable.tbl
##  (downloaded from http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22, columns: Entry; Entry name; Gene names(primary))


bash Code/preparePathStructures.sh

#######################################################
## 1.) DHS download and processing
#######################################################
#######################################################
## Code/DataFetch/getDHSfdr0.01data.sh
## input:
## output:Data/DHSfdr0.01/wgEncode(Uw|Duke)Dnase(\w+).fdr01peaks.hg19.bb.bed
##
## explanation: Downloads DHS data
#######################################################
bash Code/DataFetch/getDHSfdr0.01data.sh 

#######################################################
## Code/DHSprocessing/averageQuantileNormalizeDHS.R
##
## input: Data/DHSfdr0.01/$1.bed
## output: interimData/DHSnormalized/Normalized.$1.bed
##
## explanation: cuts all DHS-files such that 'cutoff'-top
## peaks are kept. Then quantile Normalizes signal
## where a given quantile is the mean of the quantiles of
## all underlying empirical distributions 
## (only the 'cutoff'-top feature is later used, the script is
## here for legacy purposes.)
#######################################################
Rscript  Code/DHSprocessing/averageQuantileNormalizeDHS.R

#######################################################
## 2.) Expression download and processing
#######################################################
#######################################################
##IMPORTANT: download CEL files for GEO series GSE19090 and GSE15805 manually and put into
## Data/CELS/.+ 
#######################################################
#######################################################
## Code/getIntersectionOfFiles.R
## input: Data/CELs/.+ 
## input: Data/DHSfdr0.01/wgEncode(Uw|Duke)Dnase(\w+).fdr01peaks.hg19.bb.bed
## output: interimData/fileNameTables.RDat
##
## explanation: creates a cell line mapping file to annotate the experiments
## to their respective cell line.  
#######################################################
Rscript Code/getIntersectionOfFiles.R

#######################################################
## Code/ExpressionProcessing/makeExpressionMatrixDirect.R 
## input: interimData/OverlappedCELs/.+
## output: interimData/ram_my_exprs.RDat
##
## explanation: aggregates CEL files and 
## uses rma-package: http://bmbolstad.com/misc/ComputeRMAFAQ/ComputeRMAFAQ.html
## RMA is the Robust Multichip Average. 
## It consists of three steps: a background adjustment,
## quantile normalization (see the Bolstad et al reference) and finally summarization.
## Some references (currently published) for the RMA methodology are:
## Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15 
##
### only performs analysis for CEL-files that have a corresponding dhs-file.
#######################################################
Rscript Code/ExpressionProcessing/makeExpressionMatrixDirect.R

#######################################################
## Code/ExpressionProcessing/prepareProbesetToGeneSymbolMappingDirect.R
##
## input: Data/mapping_71.txt
## input: Data/martEnsg2hgnc.txt
## output:interimData/probesetGenenameMapDirect.tbl
##
## explanation: maps the core-probeset-ids  to gene symbols.
## by intersecting mapping refseq-exon annotation to AffyExonProbesetCore Annotation.
#######################################################
Rscript Code/ExpressionProcessing/prepareProbesetToGeneSymbolMappingDirect.R

#######################################################
## Code/ExpressionProcessing/mapGeneNamesOntoExpressionsDirect.R
##
## input: interimData/expressionAsRamDirect.RDat
## input: interimData/probesetGenenameMapDirect.tbl
## output: interimData/expressionUnaveragedDirect.RDat
##
## explanation: maps GeneNames onto probesets.
## keeps only 1-to-1 mappings.
#######################################################
Rscript Code/ExpressionProcessing/mapGeneNamesOntoExpressionsDirect.R

#######################################################
## Code/ExpressionProcessing/averageExpressionDirect.R
##
## input: interimData/expressionUnaveragedDirect.RDat
## input: interimData/fileNameTables.RDat
## output: interimData/overallAveragedExpressionDirect.RDat
##
## explanation: average expression values across plates.
#######################################################
Rscript Code/ExpressionProcessing/averageExpressionDirect.R

#######################################################
## 3.) Motif download and processing
#######################################################
#######################################################
## Code/DataFetch/downloadMotifData.sh
## input:
## output:Data/motifData/HOCOMOCOv9_AD.meme
##
## explanation: downloads  motif files. 
#######################################################
bash Code/DataFetch/downloadMotifData.sh

#######################################################
## Code/DataFetch/downloadMotifData.sh
## input:
## output:Data/genomicSequence/chr(\d+).fa
##
## explanation: downloads genomic sequence.
#######################################################
bash Code/DataFetch/getGenomeData.sh

#######################################################
## Code/motifProcessing/call_fimo_HOCOMOCO.sh
##
## uses: Code/motifProcessing/run_fimo.sh
## input: Data/genomicSequence/chr(\d+).fa 
## input: Data/motifData/HOCOMOCOv9_AD.meme
## output: interimData/motifInstances/all_motifTsHOCOMOCO_(\d+)
##
## explanations: runs fimo over HOCOMOCO motifs.
#######################################################
bash  Code/motifProcessing/call_fimo_HOCOMOCO.sh

#######################################################
## Code/motifProcessing/splitMotifsHOCOMOCO.sh
##
## uses: interimData/motifInstances/HOCOMOCOSPLITTED/
## input: interimData/motifInstances/all_motifTsHOCOMOCO_(\d+)
## output: interimData/motifInstances/HOCOMOCOBEDS/CEBPD_f1.bed (for example)
## explanations: splits large fimo output files.
#######################################################
Code/motifProcessing/splitMotifsHOCOMOCO.sh

#######################################################
## Code/motifProcessing/prepareHighConfidenceMotifs.sh
##
## input:  interimData/motifInstances/HOCOMOCOBEDS/CEBPD_f1.bed (for example)
## output: interimData/motifInstances/HOCOMOCOHIGHCONF/CEBPD_f1_HighConf.bed (for example)
## explanations: keeps only motifs above 1E-5 (exist only for legacy reasons).
#######################################################
bash Code/motifProcessing/prepareHighConfidenceMotifs.sh

#######################################################
## Code/motifProcessing/make_motifDHS_matrixDirect.R
##
## input: interimData/motifInstances/HOCOMOCOHIGHCONF/CEBPD_f1_HighConf.bed (for example)
## input: interimData/DHSnormalized/Normalized.wgEncode(Uw|Duke)Dnase(\w+).fdr01peaks.hg19.bb.bed
## input: interimData/fileNameTables.RDat
## output: interimData/unnormedMotifActivityDirect.RDat
## explanations: make raw motif activity matrix.
#######################################################
Rscript Code/motifProcessing/make_motifDHS_matrixDirect.R

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
bash Code/motifProcessing/prepareMotifClusters.sh

#######################################################
## Code/mapResultsToTFClusters.R
##
## input: interimData/unnormedMotifActivityDirect.RDat
## input: Data/UniProtGeneSymbolTable.tbl
## input: interimData/tfClassifications.txt
##
## output: interimData/motifFamilyTable.RData
## output: interimData/completeMotifFamilyTable.RData
 
## explanation: maps motif-ids to geneIds via  uniprotIds
## (gives almost complete coverage). additionally,
## gives for each factor, to which subfamily and family it belongs.
## completeMotifFamilytTable.RData contains information about
## all tfs given in wingender; not just the ones that we have
## a motif for. 
#######################################################
Rscript Code/mapResultsToTFClusters.R

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
## input: interimData/allRes1.RDat
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
## input: interimData/alldf1SubFamWithVariance.RDat
## output: PaperDocs/Images/pioneerEnrichment.svg (supFig 5)
#######################################################
Rscript Code/makeFigures/addVarianceEffectToSubFamilyTable.R

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
## input: interimData/alldf1SubFamWithVariance.RDat
## output: PaperDocs/Images/showVarianceEffectSubfamily.svg
## output: aperDocs/supplementaryTable1.txt
## (supFig 3)
######################################################
Rscript Code/makeFigures/makeVarianceEffectPlot.R

#######################################################
#######################################################
### ChIP-Seq analysis
## explanation: prepare plot of Fig3
#######################################################
#######################################################
## Code/chipSeq/downloadData.sh
##
## output: Data/chipSeq/Haib/wgEncodeHaib(.+).broadPeak
## output: Data/chipSeq/Sydh/wgEncodeSydh(.+).broadPeak
## output: Data/chipSeq/wgEncodeHaibTfbsPrepared.txt
## output: Data/chipSeq/wgEncodeSydhTfbsPrepared.txt
#######################################################
Rscript Code/chipSeq/downloadData.sh

####################################################a###
## Code/chipSeq/prepareEnrichmentMatricesFinal.R
##
## uses Code/chipSeq/processNames.R
## input: interimData/fileNameTables.RDat
## output: interimData/motifInstances/HOCOMOCOHIGHCONF/CEBPD_f1_HighConf.bed (for example)
## output: Data/chipSeq/Haib/wgEncodeHaib(.+).broadPeak
## output: Data/chipSeq/Sydh/wgEncodeSydh(.+).broadPeak
#######################################################
Rscript Code/chipSeq/prepareEnrichmentMatricesFinal.R

####################################################a###
## Code/chipSeq/prepareEnrichmentMatricesSubfamilyFinal.R
##
## uses Code/chipSeq/processNames.R
## input: interimData/fileNameTables.RDat
## input: interimData/motifFamilyTable.RData
## input: interimData/DHSnormalized/Normalized.wg.$1.bed
## input: interimData/unnormedMotifActivityDirect.RDat
## input: interimData/motifInstances/HOCOMOCOHIGHCONF/CEBPD_f1_HighConf.bed (for example)
## output: interimData/chipSeqEnrichment.RDat
#######################################################
Rscript Code/chipSeq/prepareEnrichmentMatricesSubfamilyFinal.R

####################################################a###
## Code/chipSeq/prepareEnrichmentMatricesSubfamilyFinal.R
##
## uses Code/chipSeq/processNames.R
## input: interimData/fileNameTables.RDat
## input: interimData/motifFamilyTable.RData
## input: interimData/DHSnormalized/Normalized.wg.$1.bed
## input: interimData/unnormedMotifActivityDirect.RDat
## input: interimData/motifInstances/HOCOMOCOHIGHCONF/CEBPD_f1_HighConf.bed (for example)
## output: interimData/chipSeqEnrichmentSubfamily.RDat
#######################################################
Rscript Code/chipSeq/prepareEnrichmentMatricesSubfamilyFinal.R

########################################################
## Code/chipSeq/makeBoxPlotsFinal.R
##
## input: interimData/chipSeqEnrichment.RDat
## input: interimData/chipSeqEnrichmentSubfamily.RDat
## output:PaperDocs/Images/chipSeqEnrichment.svg (Fig 3)
#######################################################
Rscript Code/chipSeq/makeBoxPlotsFinal.R
#######################################################
### END: ChIP-Seq analysis
#######################################################
#######################################################
### PIQ-results comparison
#######################################################
#######################################################
## Code/compareToPIQ/prepareCombineCSRandPIQplot.R
##
## uses Code/makeFigures/helperFunctionsFinal.R
## input: Data/piqResultsTable.txt
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/alldf1SubFam.RDat
## output: PaperDocs/Images/combineCSRandPIQ.pdf
#######################################################
Rscript Code/compareToPIQ/prepareCombineCSRandPIQplot.R
#######################################################
#######################################################
### END: PIQ-results comparsion
#######################################################
