prepareNormalizedMatrices=function(nrOfPCsMotif,minimalOverlapCount){
#######################
    library(data.table)
    library(ggplot2)
    library(reshape2)
    source("Code/normalizationFunctions.R")
    load("interimData/unnormedH3K4MotifActivityDirect.RDat")
#######################
    ## remove motifs below an average overlap count
    unnormedMotifActivity=unnormedMotifActivity[rowMeans(unnormedMotifActivity)>minimalOverlapCount,]
#######################
    load("interimData/overallAveragedExpressionDirect.RDat")
    activityCellNames=colnames(unnormedMotifActivity)
    expressionCellNames=colnames(averagedExpression)
    commonNames=intersect(activityCellNames,expressionCellNames)
    expressionMatch=match(commonNames,expressionCellNames)
    activityMatch=match(commonNames,activityCellNames)

    motifActivity=unnormedMotifActivity[,activityMatch]
    unnormedExpression=averagedExpression[,expressionMatch]

    if(sum(colnames(motifActivity)!=colnames(unnormedExpression))>0){
        print("mistake in column order")
    }

    scaledExpression=scalerRound(unnormedExpression,5)
    normedExpression=scaledExpression

    motifActivity=scalerRound(motifActivity,5)
    if (nrOfPCsMotif>0){   
        topVects=getFirstEig(motifActivity,nrOfPCsMotif)
        motifActivity=removeVects(motifActivity,topVects)
        topVectsMotif=topVects
    }
    motifActivityNotQuantileNormalized=motifActivity
    motifActivity=QQnormalizeMatRow(motifActivity)
#############################
    bothMats=list(motifActivity,normedExpression,motifActivityNotQuantileNormalized,unnormedExpression)
    saveStr=paste("interimData/bothMatDirectsH3K4DeterministicPCMotif",nrOfPCsMotif,"PCsRemoved.RDat",sep="")
    save(bothMats,file=saveStr)
}


prepareNormalizedMatrices(1,100)
prepareNormalizedMatrices(0,100)
