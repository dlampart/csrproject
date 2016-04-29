library(data.table)

getZaretPioneerSubFamilies=function(){
    out=c("3.3.1.1","3.1.10.5","4.1.1.2","2.3.1.2","2.2.1.1","1.2.2.2","3.2.1.1","3.5.2.5","6.3.1.0","1.2.5.2")
    return(out)
}

getMaximalPossibleRankForSecond=function(){
    load("interimData/allRes1.RDat")
    load("interimData/completeMotifFamilyTable.RData")
    motifFamilyTableWithExpr=motifFamilyTable[is.element(geneSymbol, rownames(allRes1[[1]])),]
    nrOfSubFamMembers=motifFamilyTableWithExpr[,length(geneSymbol),by=subFamily]
    len=sum(nrOfSubFamMembers[,V1>1])
    return(len)
}

getMaximalPossibleRank=function(){
    load("interimData/allRes1.RDat")
    load("interimData/completeMotifFamilyTable.RData")
    len=length(unique(motifFamilyTable[is.element(geneSymbol,rownames(allRes1[[1]])),subFamily]))
    return(len)
}

getRank=function(res,resName,motifFamilyTable){
        fun=function(E,A){
            return(E[which.min(A)])
        }
        topName=rownames(res)[which.max(res$chi_sq)][1]
        currentSubFamily=motifFamilyTable[`UniProtKB-ID`==resName,subFamily]    
        pvals=1-pchisq(res$chi_sq,1)
        pvalRes=data.table(pvals=pvals,geneSymbol=rownames(res))
        setkey(pvalRes,geneSymbol)
        setkey(motifFamilyTable,geneSymbol)
        fused=merge(motifFamilyTable,pvalRes)
        WW=fused[,list(length(geneSymbol),min(pvals)),by=subFamily]
        WW2=WW[,list(subFamily,rank(V2*V1))]
        currentRank=WW2[subFamily==currentSubFamily,V2]
        if(length(currentRank)!=1){
            return(NULL)
        }
        allBestGeneNames=fused[,fun(geneSymbol,pvals),subFamily]
        return(data.table(rank=currentRank,topGene=allBestGeneNames[subFamily==currentSubFamily,V1],topGeneOverall=topName))
}


makeIndexListFromNamesList=function(bothMats,nameSetList){
#######################################################
## running fast_lmm_reg max and sum save
## input: bothMats (input matrices)
## input: nameSetList (list of name Sets
## output: indexSetList (list of indices to run  fast_lmm_reg (max and sum))
#######################################################
    geneNames=rownames(bothMats[[2]])
    indexSetList=list()
    for(i in c(1:length(nameSetList))){
        inds=match(nameSetList[[i]],geneNames)
        indexSetList[[i]]=inds[!is.na(inds)]
    }
    names(indexSetList)=names(nameSetList)
    return(indexSetList)
}

runFastLmmRegMaxAndSum=function(bothMats,indexList){
#######################################################
## running fast_lmm_reg max and sum save
## input: bothMats (input matrices)
## input: indexList (list of indices to run  fast_lmm_reg (max and sum))
## output: allResults : list of results.    
#######################################################
    library(multicore)
    source("Code/fast_lmm_group_reg.R")
    tfNames=rownames(bothMats[[1]])
    len=length(tfNames)
    motifActivityList=list()
    for(i in c(1:len)){
        motifActivityList[[tfNames[i]]]=bothMats[[1]][i,]
    }    
    run_fast_lmm_group_reg=function(motifActivity){
        wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,regFactor=Inf)
        return(wer$pval)
    }
    run_fast_lmm_group_regMax=function(motifActivity){
        wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,regFactor=Inf,useMax=T)
        return(wer$pval)
    }
    resultListSum=mclapply(motifActivityList,run_fast_lmm_group_reg,mc.cores=25)
    resultListMax=mclapply(motifActivityList,run_fast_lmm_group_regMax,mc.cores=25)
    makeMat=function(resultList){
        resultsMat=matrix(NA,nrow=length(resultList),ncol=length(resultList[[1]]))
        tfNames=rownames(bothMats[[1]])
        len=length(resultsMat[,1])
        for(i in c(1:len)){
            resultsMat[i,]=resultList[[i]]
        }
        rownames(resultsMat)=tfNames
        colnames(resultsMat)=names(indexList)
        return(resultsMat)
    }
    allResults=list()
    allResults[["resultMatSum"]]=makeMat(resultListSum)
    allResults[["resultMatMax"]]=makeMat(resultListMax)
    return(allResults)
}



runFastLmmRegMax=function(bothMats,indexList, covariates=NULL,myK=NULL){
#######################################################
## running fast_lmm_reg max and sum save
## input: bothMats (input matrices)
## input: indexList (list of indices to run  fast_lmm_reg (max))
## output: allResults : list of results.    
#######################################################
    library(multicore)
    source("Code/fast_lmm_group_reg.R")
    source("Code/fast_lmm_covariates.R")
    tfNames=rownames(bothMats[[1]])
    len=length(tfNames)
    motifActivityList=list()
    for(i in c(1:len)){
        motifActivityList[[tfNames[i]]]=bothMats[[1]][i,]
    }    
    run_fast_lmm_group_reg=function(motifActivity){
        wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,regFactor=Inf)
        return(wer$pval)
    }
    if(is.null(covariates)){
        run_fast_lmm_group_regMax=function(motifActivity){
            wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,regFactor=Inf,useMax=T)
            return(wer$pval)
        }    
    }else{
        run_fast_lmm_group_regMax=function(motifActivity){

                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=FALSE,useMax=T,covariateMat=covariates)
 #                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=TRUE,covariateMat=covariates)
            return(wer$pval)
        }
    }
#browser()
                                        #    resultListSum=mclapply(motifActivityList,run_fast_lmm_group_reg,mc.cores=15)
   resultListMax=mclapply(motifActivityList,run_fast_lmm_group_regMax,mc.cores=25)
#        resultListMax=lapply(motifActivityList,run_fast_lmm_group_regMax)
    makeMat=function(resultList){
        resultsMat=matrix(NA,nrow=length(resultList),ncol=length(resultList[[1]]))
        tfNames=rownames(bothMats[[1]])
        len=length(resultsMat[,1])
        for(i in c(1:len)){
            resultsMat[i,]=resultList[[i]]
        }
        rownames(resultsMat)=tfNames
        colnames(resultsMat)=names(indexList)
        return(resultsMat)
    }
    allResults=list()
 #   allResults[["resultMatSum"]]=makeMat(resultListSum)
    allResults[["resultMatMax"]]=makeMat(resultListMax)
    return(allResults)
}

runFastLmmRegMaxAllResults=function(bothMats,indexList, covariates=NULL,myK=NULL){
#######################################################
## running fast_lmm_reg max and sum save
## input: bothMats (input matrices)
## input: indexList (list of indices to run  fast_lmm_reg (max))
## output: allResults : list of results.    
#######################################################
    library(multicore)
    source("Code/fast_lmm_group_reg.R")
    source("Code/fast_lmm_covariates.R")
    tfNames=rownames(bothMats[[1]])
    len=length(tfNames)
    motifActivityList=list()
    for(i in c(1:len)){
        motifActivityList[[tfNames[i]]]=bothMats[[1]][i,]
    }    
    run_fast_lmm_group_reg=function(motifActivity){        
        wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=myK,indexList=indexList,regFactor=Inf)
        return(wer$pval)
    }
    if(is.null(covariates)){
        run_fast_lmm_group_regMax=function(motifActivity){
            wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=myK,indexList=indexList,regFactor=Inf,useMax=T,renorm=F)            
            return(wer$pval)
        }    
    }else{
        run_fast_lmm_group_regMax=function(motifActivity){

                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=FALSE,useMax=T,covariateMat=covariates)
 #                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=TRUE,covariateMat=covariates)
            return(wer$pval)
        }
    }

        getOnlyDeltasWrapper=function(motifActivity){
                       wer=getOnlyDelta(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,regFactor=Inf,useMax=T)
 #                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=TRUE,covariateMat=covariates)
            return(wer)
        }

                                        #browser()
                                        #    resultListSum=mclapply(motifActivityList,run_fast_lmm_group_reg,mc.cores=15)

#    resultListMax=lapply(motifActivityList,run_fast_lmm_group_regMax)
   resultListMax=mclapply(motifActivityList,run_fast_lmm_group_regMax,mc.cores=25)

    deltas=mclapply(motifActivityList,getOnlyDeltasWrapper,mc.cores=25)
#        deltas=lapply(motifActivityList,getOnlyDeltasWrapper)
#        browser()
#        resultListMax=lapply(motifActivityList,run_fast_lmm_group_regMax)
    makeMat=function(resultList){
        resultsMat=matrix(NA,nrow=length(resultList),ncol=length(resultList[[1]]))
        tfNames=rownames(bothMats[[1]])
        len=length(resultsMat[,1])
        for(i in c(1:len)){
            resultsMat[i,]=resultList[[i]]
        }
        rownames(resultsMat)=tfNames
        colnames(resultsMat)=names(indexList)
        return(resultsMat)
    }
    allResults=list()
 #   allResults[["resultMatSum"]]=makeMat(resultListSum)
    allResults[["resultMatMax"]]=makeMat(resultListMax)
    allResults[["deltas"]]=deltas
    return(allResults)
}

runFastLmmRegGroupAllResults=function(bothMats,indexList, covariates=NULL,myK=NULL){
#######################################################
## running fast_lmm_reg max and sum save
## input: bothMats (input matrices)
## input: indexList (list of indices to run  fast_lmm_reg (max))
## output: allResults : list of results.    
#######################################################
    library(multicore)
    source("Code/fast_lmm_group.R")
    source("Code/fast_lmm_covariates.R")
    tfNames=rownames(bothMats[[1]])
    len=length(tfNames)
    motifActivityList=list()
    for(i in c(1:len)){
        motifActivityList[[tfNames[i]]]=bothMats[[1]][i,]
    }    
    run_fast_lmm_group=function(motifActivity){        
        wer=fast_lmm_group(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=myK,indexList=indexList)
        return(wer$pval)
    }
    if(is.null(covariates)){
        run_fast_lmm_group_regMax=function(motifActivity){
            wer=fast_lmm_group(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=myK,indexList=indexList)            
            return(wer$pval)
        }    
    }else{
        run_fast_lmm_group_regMax=function(motifActivity){

                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=FALSE,useMax=T,covariateMat=covariates)
 #                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=TRUE,covariateMat=covariates)
            return(wer$pval)
        }
    }

        getOnlyDeltasWrapper=function(motifActivity){
                       wer=getOnlyDelta(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,regFactor=Inf,useMax=T)
 #                       wer=fast_lmm_covariates(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,fullLikelihood=TRUE,covariateMat=covariates)
            return(wer)
        }

                                        #browser()
                                        #    resultListSum=mclapply(motifActivityList,run_fast_lmm_group_reg,mc.cores=15)

#    browser()
#    resultListMax=lapply(motifActivityList,run_fast_lmm_group_regMax)
   resultListMax=mclapply(motifActivityList,run_fast_lmm_group_regMax,mc.cores=25)

    deltas=mclapply(motifActivityList,getOnlyDeltasWrapper,mc.cores=25)
#        deltas=lapply(motifActivityList,getOnlyDeltasWrapper)
#        browser()
#        resultListMax=lapply(motifActivityList,run_fast_lmm_group_regMax)
    makeMat=function(resultList){
        resultsMat=matrix(NA,nrow=length(resultList),ncol=length(resultList[[1]]))
        tfNames=rownames(bothMats[[1]])
        len=length(resultsMat[,1])
        for(i in c(1:len)){
            resultsMat[i,]=resultList[[i]]
        }
        rownames(resultsMat)=tfNames
        colnames(resultsMat)=names(indexList)
        return(resultsMat)
    }
    allResults=list()
 #   allResults[["resultMatSum"]]=makeMat(resultListSum)
    allResults[["resultMatMax"]]=makeMat(resultListMax)
    allResults[["deltas"]]=deltas
    return(allResults)
}


runFastLmmRegMaxWithSecMat=function(bothMats,indexList, covariates=NULL){
#######################################################
## running fast_lmm_reg max and sum save
## input: bothMats (input matrices)
## input: indexList (list of indices to run  fast_lmm_reg (max))
## output: allResults : list of results.    
#######################################################
    library(multicore)
#    source("Code/fast_lmm_group_reg.R")
    source("Code/fast_lmmEM.R")
    tfNames=rownames(bothMats[[1]])
    len=length(tfNames)
    motifActivityList=list()
    for(i in c(1:len)){
        motifActivityList[[tfNames[i]]]=bothMats[[1]][i,]
    }
        run_fast_lmmEM=function(motifActivity){
            wer=fast_lmmEM(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=NULL,indexList=indexList,useMax=T,my_K2=cov(bothMats[[1]])) 
#                          wer=fast_lmm_group_reg(my_y=as.matrix(motifActivity), my_X=t(bothMats[[2]]), my_K=cov(bothMats[[1]]),indexList=indexList,useMax=T)
            return(wer$pval)
        }
    
  resultListMax=mclapply(motifActivityList,run_fast_lmmEM,mc.cores=25)
#resultListMax=lapply(motifActivityList,run_fast_lmmEM)
    makeMat=function(resultList){
        resultsMat=matrix(NA,nrow=length(resultList),ncol=length(resultList[[1]]))
        tfNames=rownames(bothMats[[1]])
        len=length(resultsMat[,1])
        for(i in c(1:len)){
            resultsMat[i,]=resultList[[i]]
        }
        rownames(resultsMat)=tfNames
        colnames(resultsMat)=names(indexList)
        return(resultsMat)
    }
    allResults=list()
 #   allResults[["resultMatSum"]]=makeMat(resultListSum)
    allResults[["resultMatMax"]]=makeMat(resultListMax)
    return(allResults)
}



makeSetListWingender=function(motifFamilyTable,bothMats,quotedStr,returnNames=F){
#######################################################
## input: motifFamilyTable (from Code/mapResultsToTFClusters.R)
## input:  bothMats (from Code/prepareNormalizedMatrices.R)
## input:  quotedStr: (either quote(subFamily) or quote(family)
## output: indexSetList: list of indices
#######################################################
    geneNames=rownames(bothMats[[2]])
    motifFamilyTable=as.data.table(motifFamilyTable)
    indexSetListNames=list()
    indexSetList=list()
    allSubFamilies=motifFamilyTable[,unique(eval(quotedStr))]
    len=length(allSubFamilies)
    for(i in c(1:len)){
        currentGeneNames=motifFamilyTable[allSubFamilies[i]==eval(quotedStr),geneSymbol]
        currentGeneNames=intersect(currentGeneNames,geneNames)
        if(length(currentGeneNames)>0){
            indexSetListNames[[allSubFamilies[i]]]=currentGeneNames
            indexSetList[[allSubFamilies[i]]]=which(is.element(geneNames,currentGeneNames))
        }    
    }
    if(returnNames){
        return(indexSetListNames)
    }
    return(indexSetList)
}





repListNamesToMatchElNr=function(myList){
    ss=unlist(lapply(myList,length))
    myNames=names(ss)
    ret=rep(myNames,ss)
    return(ret)
}

findSetOfMotif=function(bothMats,indexSetList,motifFamilyTable){
#######################################################
## input: motifFamilyTable (from Code/mapResultsToTFClusters.R)
## input:  bothMats (from Code/prepareNormalizedMatrices.R)
## input: indexSetList: list of indices (index refer to rownames(bothMats[[2]])
## output: setNames: which set name does the motif Belong to
#######################################################
    geneNames=rownames(bothMats[[2]])
    tfNames=rownames(bothMats[[1]])
    len=length(tfNames)
    tfProteins=sub("_..$","",tfNames)
    motifInSubFamilyIndex=c(NA,len)
    tfGenes=motifFamilyTable[match(tfProteins,`UniProtKB-ID`),geneSymbol]
    myNames=repListNamesToMatchElNr(indexSetList)
    inds=unlist(indexSetList,use.names=F)
    tfGeneNamesMatched=geneNames[inds]
    myMatchInds=match(tfGenes,tfGeneNamesMatched)
    setNames=myNames[myMatchInds]
    return(setNames)
}

findSetOfMotifNovel=function(bothMats,indexSetList,motifFamilyTable,novelMotifsMapped,pvalCutoff){
#######################################################
## input: motifFamilyTable (from Code/mapResultsToTFClusters.R)
## input:  bothMats (from Code/prepareNormalizedMatrices.R)
## input: indexSetList: list of indices (index refer to rownames(bothMats[[2]])
## input: novelMotifsMapped: table of mapping between Novel motifs and known motifs (with TOMTOM-pvalues)
## input: pvalCutoff: param gives which mappings in novelMotifsMapped  are to be considered.
## output: setNames: which set name does the motif Belong to
#######################################################
    novelMotifsMapped=novelMotifsMapped[pvalue<pvalCutoff,]
    geneNames=rownames(bothMats[[2]])
    tfNames=sub(".bed","",rownames(bothMats[[1]]))
    len=length(tfNames)
    tfGenes=novelMotifsMapped[match(tfNames,`#Query ID`),geneSymbols]
    motifInSubFamilyIndex=c(NA,len)
#    tfGenes=motifFamilyTable[match(tfProteins,`UniProtKB-ID`),geneSymbol]
    myNames=repListNamesToMatchElNr(indexSetList)
    inds=unlist(indexSetList,use.names=F)
    tfGeneNamesMatched=geneNames[inds]
    myMatchInds=match(tfGenes,tfGeneNamesMatched)
    setNames=myNames[myMatchInds]
    return(setNames)
}

plotEnrichment=function(plotFileName,resultsMat,motifInSetIndex,plotTitle=""){   
#######################################################
## running plotting scaled rank histogram (rank is calculated as the rank of set containing correct motif)
## input: plotFileName: name of output plot
## input: resultsMat (submatrix derived from runFastLmmRegMaxAndSum)
## input: motifInSetIndex (index of set containing the motif (NA if not found)
## output: rank vector
#######################################################
    remove=is.na(motifInSetIndex)
    resultsMat=resultsMat[!remove,]
    motifInSetIndex=motifInSetIndex[!remove]    
    resultsMatRanked=apply(resultsMat,1,rank)
    resultsMatUnif=resultsMatRanked/length(resultsMat[1,])
    len=length(motifInSetIndex)
    rankOfMotif=rep(NA,len)
    for(i in c(1:len)){
        ind=motifInSetIndex[i]
        if (length(ind)>0)
            rankOfMotif[i]=resultsMatUnif[ind,i]
    }
    data=data.table(rankOfMotif,rankOfMotif)
    ggplot(data=data,aes(x=rankOfMotif))+geom_histogram()+xlim(c(0,1))+xlab("theoretical")+ylab("empirical")+ggtitle(plotTitle)+theme(axis.text=element_text(size=13),axis.text=element_text(size=16,face="bold"))
    ggsave(plotFileName)
    return(rankOfMotif*length(resultsMat[1,]))
}

plotEnrichmentAndSave=function(plotFileName,resultsMat,motifInSetIndex,plotTitle=""){   
#######################################################
## running plotting scaled rank histogram (rank is calculated as the rank of set containing correct motif)
## input: plotFileName: name of output plot
## input: resultsMat (submatrix derived from runFastLmmRegMaxAndSum)
## input: motifInSetIndex (index of set containing the motif (NA if not found)
## output: rank vector
#######################################################
    remove=is.na(motifInSetoIndex)
    resultsMat=resultsMat[!remove,]
    motifInSetIndex=motifInSetIndex[!remove]    
    resultsMatRanked=apply(resultsMat,1,rank)
    resultsMatUnif=resultsMatRanked/length(resultsMat[1,])
    len=length(motifInSetIndex)
    rankOfMotif=rep(NA,len)
    for(i in c(1:len)){
        ind=motifInSetIndex[i]
        if (length(ind)>0)
            rankOfMotif[i]=resultsMatUnif[ind,i]
    }
    data=data.table(rankOfMotif=rankOfMotif,names=rownames(resultsMat))
    ggplot(data=data,aes(x=rankOfMotif))+geom_histogram()+xlim(c(0,1))+xlab("theoretical")+ylab("empirical")+ggtitle(plotTitle)+theme(axis.text=element_text(size=13),axis.text=element_text(size=16,face="bold"))
    ggsave(plotFileName)
    browser()
    return(rankOfMotif*length(resultsMat[1,]))
}


plotEnrichmentPaneled=function(plotFileName,resultsMat,motifInSetIndex,superclassNames){   
#######################################################
## running plotting scaled rank histogram (rank is calculated as the rank of set containing correct motif)
## input: plotFileName: name of output plot
## input: resultsMat (submatrix derived from runFastLmmRegMaxAndSum)
## input: motifInSetIndex (index of set containing the motif (NA if not found)
## output: rank vector
#######################################################

    load("interimData/unnormedMotifActivityDirect.RDat")

    remove=is.na(motifInSetIndex)
    resultsMat=resultsMat[!remove,]
    motifInSetIndex=motifInSetIndex[!remove]
    superclassNames=superclassNames[!remove]
    
    resultsMatRanked=apply(resultsMat,1,rank)
    resultsMatUnif=resultsMatRanked/length(resultsMat[1,])
    len=length(motifInSetIndex)
    rankOfMotif=rep(NA,len)
    for(i in c(1:len)){
        ind=motifInSetIndex[i]
        if (length(ind)>0)
            rankOfMotif[i]=resultsMatUnif[ind,i]
    }
    data=data.table(rankOfMotif=rankOfMotif,superclassNames=superclassNames)
    ggplot(data=data[],aes(x=rankOfMotif))+geom_histogram()+xlim(c(0,1))+facet_grid(.~superclassNames)
    ggsave(plotFileName)
    return(data)
}


makeQQplot=function(plotFileName,resultsMat,plotTitle=""){   
#######################################################
## make qqplot over all p-values:
## input: plotFileName: name of output plot
## input: resultsMat (submatrix derived from runFastLmmRegMaxAndSum)
#######################################################
    ss=sort(resultsMat)
    data=data.frame(pvals=ss,theor=c(1:length(ss))/(1+length(ss)))
    ggplot(data=data,aes(y=-log10(pvals),x=-log10(theor)))+geom_point()+geom_abline()+xlab("theoretical")+ylab("empirical")+ggtitle(plotTitle)+theme(axis.title=element_text(size=18,face="bold") ,axis.text=element_text(size=14))
    ggsave(plotFileName)
}


flagAnnotatedFamilies=function(familyCoverage,quotedStr,setList){
#######################################################
## flag sets in setlist that are given in familyCoverage (derived in Code/findAnnotatedFamilies.R)    
## input: quotedStr: flags whether family or subfamily results should be returned.
## input: setList: names of sets.
## ouput: indicator vect.    
#######################################################
    if (quotedStr==quote(family)){
        coveredFamilyIds=familyCoverage[["coveredFamilies"]]
    }
    if (quotedStr==quote(subFamily)){
        coveredFamilyIds=familyCoverage[["coveredSubFamilies"]]
    }
    return(is.element(setList,coveredFamilyIds))
}



findFDR=function(pValMatrix,fdrLevel){
#######################################################
#######################################################
    nrOfTests=length(pValMatrix[,1])*length(pValMatrix[1,])
    ee=sort(pValMatrix)
    theor=c(1:nrOfTests)/(1+nrOfTests)
    belowFDR=ee/theor<fdrLevel
    topV=max(which(belowFDR))
    return(ee[topV])        
}
