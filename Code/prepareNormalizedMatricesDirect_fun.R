prepareNormalizedMatrices=function(nrOfPCsMotif,minimalOverlapCount,inputstrMotif="interimData/unnormedMotifActivityDirect.RDat",inputstrExpr="interimData/overallAveragedExpressionDirect.RDat",savestrinfix="",normExp=TRUE){
#######################
    library(data.table)
    library(ggplot2)
    library(reshape2)
    source("Code/normalizationFunctions.R")
    load(inputstrMotif)
#######################
    ## remove motifs below an average overlap count
    unnormedMotifActivity=unnormedMotifActivity[rowMeans(unnormedMotifActivity)>minimalOverlapCount,]
#######################
    load(inputstrExpr)
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
    if(normExp){
        scaledExpression=scalerRound(unnormedExpression,5)
        normedExpression=scaledExpression
    }else{
        normedExpression=unnormedExpression
    }
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
    saveStr=paste("interimData/bothMatDirectsDeterministicPCMotif",nrOfPCsMotif,"PCsRemoved",savestrinfix,".RDat",sep="")
    save(bothMats,file=saveStr)
    print("save to")
    print(saveStr)
}
