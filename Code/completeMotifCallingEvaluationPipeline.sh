#!/bin/bash

#######################################################
## Code/motifCallingPipelineEvaluation/call_fimo_HOCOMOCO_low_conf.sh
##
## uses: Code/motifProcessing/run_fimo.sh
## input: Data/genomicSequence/chr(\d+).fa 
## input: Data/motifData/HOCOMOCOv9_AD.meme
## output: interimData/motifInstances/all_motifTsHOCOMOCO_low_conf_(\d+)
##
## explanations: runs fimo over HOCOMOCO motifs.
#######################################################
bash  Code/motifCallingPipelineEvaluation/call_fimo_HOCOMOCO_low_conf.sh

#######################################################
## Code/motifCallingPipelineEvaluation/splitMotifsRevision.sh
##
## input: interimData/motifInstances/all_motifTsHOCOMOCO_low_conf_(\d+)
## output: interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM(6|5|45)/CEBPD_f1.sorted.bed (for example)
#######################################################
bash Code/motifCallingPipelineEvaluatio/splitMotifsRevision.sh

#######################################################
## Code/motifCallingPipelineEvaluation/make_motifDHS_matrixVariableLevels.R
##
## input: interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM(6|5|45)/CEBPD_f1.sorted.bed (for example)
## input: interimData/DHSnormalized/Normalized.wgEncode(Uw|Duke)Dnase(\w+).fdr01peaks.hg19.bb.bed
## input: interimData/fileNameTables.RDat
## output: interimData/unnormedMotifActivity_motifCutoff(4.5|5|6).RDat
#######################################################
Rscript Code/motifCallingPipelineEvaluation/make_motifDHS_matrixVariableLevels.R

#######################################################
## Code/motifCallingPipelineEvaluation/prepareNormalizedMatrices_motifCutoff.R
##
## uses Code/normalizationFunctions.R
## input: interimData/unnormedMotifActivity_motifCutoff(4.5|5|6).RDat
## input: interimData/overallAveragedExpressionDirect.RDat
## output: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff(4.5|5|6).RDat
##
## explanations: normalizes both matrices and removes principal compnent with
## largest eigenvector from the motif matrix.
#######################################################
Rscript  Code/motifCallingPipelineEvaluation/prepareNormalizedMatrices_motifCutoff.R

#######################################################
## Code/makeFigures/runAllMethodsOnlyForTFsFinal.R
##
## uses: Code/makeFigures/helperFunctionsFinal.R
## uses: Code/fast_lmm.R
## uses: Code/makeFigures/runAllMethodsOnlyForTFs_fun.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff55.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff5.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff6.RDat
##
## output: interimData/allRes_55.RDat
## output: interimData/allRes_5.RDat
## output: interimData/allRes_6.RDat
##
## explanations: runs fast_lmm only for Transcription factors
#######################################################
Rscript Code/motifCallingPipelineEvaluation/runAllMethodsOnlyForTFs_cutoff.R

#######################################################
## Code/motifCallingPipelineEvaluation/prepareRankingTables_cutoff.R
##
## uses: Code/makeFigures/helperFunctionsFinal.R
## uses: Code/makeFigures/prepareRankingTables_fun.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff55.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff5.RDat
## input:  interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff6.RDat
## input: interimData/allRes_55.RDat
## input: interimData/allRes_5.RDat
## input: interimData/allRes_6.RDat
##
## output: interimData/alldf1_55.RDat
## output: interimData/alldf1_5.RDat
## output: interimData/alldf1_6.RDat
##
## explanation: prepares tables motif rank score
#######################################################
Rscript Code/motifCallingPipelineEvaluation/prepareRankingTables_cutoff.R

#######################################################
## Code/makeFigures/prepareRankingTablesPerSubFamilyFinal.R
##
## use: Code/makeFigures/helperFunctionsFinal.R
## use: Code/makeFigures/prepareRankingTablesPerSubFamily_fun.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/allRes_55.RDat
## input: interimData/allRes_5.RDat
## input: interimData/allRes_6.RDat
## input: interimData/alldf1_55.RDat
## input: interimData/alldf1_5.RDat
## input: interimData/alldf1_6.RDat
##
## output: interimData/alldf1SubFam_55.RDat
## output: interimData/alldf1SubFam_5.RDat
## output: interimData/alldf1SubFam_6.RDat
#######################################################
Rscript Code/motifCallingPipelineEvaluation/prepareRankingTablesPerSubFamily_cutoff.R

#######################################################
## Code/makeFigures/showPowerIncreasePerSubfamilyFinal.R
##
## uses Code/makeFigures/showPowerIncreasePerSubfamily_fun.R
## input: interimData/alldfMinusSubFam.RDat
## input: interimData/allRes1.RDat
## input: interimData/alldf0SubFam.RDat
## input: interimData/alldf1SubFam.RDat
## output: PaperDocs/Images/showPowerCumulativePerSubfamily.svg
#######################################################
Rscript Code/motifCallingPipelineEvaluation/showPowerIncreasePerSubfamily_cutoff.R

#######################################################
## Code/motifCallingPipelineEvaluation/showOverlap_cutoffFinal.R
##
## uses Code/makeFigures/helperFunctionsFinal.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/alldf1SubFam_55.RDat
## input: interimData/alldf1SubFam_5.RDat
## input: interimData/alldf1SubFam_6.RDat
## output: PaperDocs/Images/cutoffOverlapPlot_All.png
## output: PaperDocs/Images/cutoffOverlapPlot_All.svg
#######################################################
Rscript Code/motifCallingPipelineEvaluation/showOverlap_cutoffFinal.R

#######################################################
## Code/motifCallingPipelineEvaluation/countMotifs.R
##
## interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/
## interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM5/
## interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM6/
#######################################################
Rscript Code/motifCallingPipelineEvaluation/countMotifs.R

#######################################################
##  Code/motifCallingPipelineEvaluation/calculateMotifCutoffFinal.R
##
## input: interimData/fileNameTables.RDat
## input: interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/
## input: Data/chipSeq/
## input: interimData/DHSnormalized
## output: interimData/chipSeqVals_0.35.txt
## output: interimData/chipSeqVals_0.5.txt
## output: interimData/chipSeqVals_0.7.txt
#######################################################
Rscript Code/motifCallingPipelineEvaluation/calculateMotifCutoffFinal.R

#######################################################
##  Code/motifCallingPipelineEvaluation/make_motifDHS_matrixDirect_ChipSeq.R
##
## input: interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/CEBPD_f1.sorted.bed (for example)
## input: interimData/DHSnormalized/Normalized.wgEncode(Uw|Duke)Dnase(\w+).fdr01peaks.hg19.bb.bed
## input: interimData/fileNameTables.RDat
## output: interimData/unnormedMotifActivity_motifCutoff(4.5|5|6).RDat
#######################################################
Rscript Code/motifCallingPipelineEvaluation/make_motifDHS_matrixDirect_ChipSeq.R

#######################################################
## Code/motifCallingPipelineEvaluation/prepareNormalizedMatrices_motifCutoff_chipSeq.R
##
## uses: Code/normalizationFunctions.R
## input: interimData/unnormedMotifActivityDirectChipSeq_0.35.RDat
## input: interimData/unnormedMotifActivityDirectChipSeq_0.5.RDat
## input: interimData/unnormedMotifActivityDirectChipSeq_0.7.RDat
## input: interimData/overallAveragedExpressionDirect.RDat
## output: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff0.35_chipseq.RDat
## output: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff0.5_chipseq.RDat
## output: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff0.7_chipseq.RDat
#######################################################
Rscript Code/motifCallingPipelineEvaluation/prepareNormalizedMatrices_motifCutoff_chipSeq.R

#######################################################
## Code/makeFigures/runAllMethodsOnlyForTFsFinal.R
##
## uses: Code/makeFigures/helperFunctionsFinal.R
## uses: Code/fast_lmm.R
## uses: Code/makeFigures/runAllMethodsOnlyForTFs_fun.R
## input:interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff0.7_chipseq.RDat
## input:interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff0.5_chipseq.RDat
## input:interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff0.35_chipseq.RDat
##
## output: interimData/allRes0.7_chipseq.RDat
## output: interimData/allRes0.5_chipseq.RDat
## output: interimData/allRes0.35_chipseq.RDat
#######################################################
Rscript Code/motifCallingPipelineEvaluation/runAllMethodsOnlyForTFs_chipseq.R

#######################################################
## Code/motifCallingPipelineEvaluation/prepareRankingTables_cutoff.R
##
## uses: Code/makeFigures/helperFunctionsFinal.R
## uses: Code/makeFigures/prepareRankingTables_fun.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff55.RDat
## input: interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff5.RDat
## input:  interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved_motifCutoff6.RDat
## input: interimData/allRes_55.RDat
## input: interimData/allRes_5.RDat
## input: interimData/allRes_6.RDat
##
## output: interimData/alldf1_0.7_chipseq.RDat
## output: interimData/alldf1_0.5_chipseq.RDat
## output: interimData/alldf1_0.35_chipseq.RDat
##
## explanation: prepares tables motif rank score
#######################################################
Rscript Code/motifCallingPipelineEvaluation/prepareRankingTables_chipseq.R

#######################################################
## Code/makeFigures/prepareRankingTablesPerSubFamilyFinal.R
##
## use: Code/makeFigures/helperFunctionsFinal.R
## use: Code/makeFigures/prepareRankingTablesPerSubFamily_fun.R
## input: interimData/completeMotifFamilyTable.RData
## input: interimData/allRes_55.RDat
## input: interimData/allRes_5.RDat
## input: interimData/allRes_6.RDat
## input: interimData/alldf1_55.RDat
## input: interimData/alldf1_5.RDat
## input: interimData/alldf1_6.RDat
##
## output: interimData/alldf1SubFam_55.RDat
## output: interimData/alldf1SubFam_5.RDat
## output: interimData/alldf1SubFam_6.RDat
#######################################################
Rscript Code/motifCallingPipelineEvaluation/prepareRankingTablesPerSubFamily_chipseq.R

#######################################################
## Code/motifCallingPipelineEvaluation/showPowerIncreasePerSubfamily_chipseq.R
##
## 
##
#######################################################
Rscript Code/motifCallingPipelineEvaluation/showPowerIncreasePerSubfamily_chipseq.R


Rscript Code/motifCallingPipelineEvaluation/showOverlap_chipseqFinal.R




