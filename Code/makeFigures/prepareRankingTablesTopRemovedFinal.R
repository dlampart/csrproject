load("interimData/completeMotifFamilyTable.RData")
library(ggplot2)
library(data.table)
source("Code/makeFigures/helperFunctionsFinal.R")

getRankOfSecond=function(res,resName,motifFamilyTable){
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
###remove top genes ####
    AA=fused[,list(geneSymbol[which.min(pvals)],min(pvals)),by=subFamily]
    fused=fused[!is.element(geneSymbol,AA[,V1]),]
###remove top genes ####
    WW=fused[,list(length(geneSymbol),min(pvals)),by=subFamily]
    WW2=WW[,list(subFamily,rank(V2*V1))]
    currentRank=WW2[subFamily==currentSubFamily,V2]
    if(length(currentRank)!=1){
        return(NULL)
    }
    allBestGeneNames=fused[,fun(geneSymbol,pvals),subFamily]    
    return(data.table(rank=currentRank,topGene=allBestGeneNames[subFamily==currentSubFamily,V1],topGeneOverall=topName))
}


makeDFTopRemoved=function(allRes,motifFamilyTable){
    protName=sub("_..$","",names(allRes))
    resNames=protName[1]
    res=allRes[[1]]    
    allRanks=rep(NA,length(protName))
    resNames=protName[1]
    alldf=NULL
    alldfFirstRankRemoved=NULL
    isEmpty=rep(F,length(protName))
    isEmpty2=rep(F,length(protName))
    for(i in c(1:length(protName))){
        print(i)
        resName=protName[i]
        res=allRes[[i]]
        df=getRank(res,resName,motifFamilyTable)
        df2=getRankOfSecond(res,resName,motifFamilyTable)
        if(is.null(df) || dim(df)[1]==0){isEmpty[i]=T}else{
            df[,index:=i]
            df[,maxchi2:=max(res$chi_sq)]
        }
        if(is.null(df2) || dim(df2)[1]==0){isEmpty2[i]=T}else{        
            df2[,index:=i]
            df2[,maxchi2:=max(res$chi_sq)]
        }
        alldf=rbind(alldf,df)
        alldfFirstRankRemoved=rbind(alldfFirstRankRemoved,df2)
    }
    alldf[,geneNr:="first"]
    alldfFirstRankRemoved[,geneNr:="second"]
    alldf=rbind(alldf,alldfFirstRankRemoved)
    alldf[,motifName:=names(allRes)[index]]
    return(alldf)
}


load("interimData/allRes1TopRemoved.RDat")
load("interimData/alldf1.RDat")
names(allRes1TopRemoved)=alldf1[geneNr=="first",motifName]
alldf1TopRemoved=makeDFTopRemoved(allRes1TopRemoved,motifFamilyTable)
save(alldf1TopRemoved,file="interimData/alldf1TopRemoved.RDat")
load("interimData/allRes1TopRemovedShuffled.RDat")
load("interimData/alldf1.RDat")
names(allRes1TopRemovedShuffled)=alldf1[geneNr=="first",motifName]
alldf1TopRemovedShuffled=makeDFTopRemoved(allRes1TopRemovedShuffled,motifFamilyTable)
save(alldf1TopRemovedShuffled,file="interimData/alldf1TopRemovedShuffled.RDat")


