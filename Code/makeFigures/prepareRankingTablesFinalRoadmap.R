load("interimData/completeMotifFamilyTable.RData")
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
###END remove top genes ####
    WW=fused[,list(length(geneSymbol),min(pvals)),by=subFamily]
    WW2=WW[,list(subFamily,rank(V2*V1))]
    currentRank=WW2[subFamily==currentSubFamily,V2]
    if(length(currentRank)!=1){
        return(NULL)
    }
    allBestGeneNames=fused[,fun(geneSymbol,pvals),subFamily]    
    return(data.table(rank=currentRank,topGene=allBestGeneNames[subFamily==currentSubFamily,V1],topGeneOverall=topName))
}

makeDF=function(bothMats,allRes,motifFamilyTable){
    protName=sub("_..$","",rownames(bothMats[[1]]))
    resNames=protName[1]
    res=allRes[[1]]    
    allRanks=rep(NA,length(protName))
    resNames=protName[1]
    alldf=NULL
    alldfFirstRankRemoved=NULL
    isEmpty=rep(F,length(protName))
    isEmpty2=rep(F,length(protName))
    for(i in c(1:length(protName))){
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
    alldf[,motifName:=rownames(bothMats[[1]])[index]]
    return(alldf)
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

removeTopGenesFromAllRes=function(allRes,motifFamilyTable,bothMats){
    allRes2=allRes
    motifFamilyAllResOrder=data.table(motifFamilyTable[match(rownames(allRes[[i]]),motifFamilyTable[,geneSymbol]),],chi_sq=allRes[[i]]$chi_sq)
    fun=function(x,y){
        return(y[which.max(x)])
    }
    for(i in c(1:length(allRes))){
        motifFamilyAllResOrder[,chi_sq:=allRes[[i]]$chi_sq]
        aa=motifFamilyAllResOrder[,fun(chi_sq,tfId),by=subFamily][,V1]
        motifFamilyAllResOrder[is.element(tfId,aa),chi_sq:=0]
        allRes2[[i]]$chi_sq=motifFamilyAllResOrder[,chi_sq]
    }
    return(allRes2)
}

#load("interimData/allResMinus.RDat")
#load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
#alldfMinus=makeDF(bothMats,allResMinus,motifFamilyTable)
#save(alldfMinus,file="interimData/alldfMinus.RDat")
#load("interimData/allRes0.RDat")
#load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
#alldf0=makeDF(bothMats,allRes0,motifFamilyTable)
#save(alldf0,file="interimData/alldf0.RDat")
#load("interimData/allRes1.RDat")
#load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
#alldf1=makeDF(bothMats,allRes1,motifFamilyTable)
#save(alldf1,file="interimData/alldf1.RDat")

#load("interimData/allRes1TopRemoved.RDat")
#load("interimData/alldf1.RDat")
#names(allRes1TopRemoved)=alldf1[geneNr=="first",motifName]
#alldf1TopRemoved=makeDFTopRemoved(allRes1TopRemoved,motifFamilyTable)
#save(alldf1TopRemoved,file="interimData/alldf1TopRemoved.RDat")
#load("interimData/allRes1TopRemovedShuffled.RDat")
#load("interimData/alldf1.RDat")
#names(allRes1TopRemovedShuffled)=alldf1[geneNr=="first",motifName]
#alldf1TopRemovedShuffled=makeDFTopRemoved(allRes1TopRemovedShuffled,motifFamilyTable)
#save(alldf1TopRemovedShuffled,file="interimData/alldf1TopRemovedShuffled.RDat")

#load("interimData/allRes1.RDat")
#load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
#alldf1=makeDF(bothMats,allRes1,motifFamilyTable)
#allResTopIgnored=removeTopGenesFromAllRes(allRes1,motifFamilyTable)
#save(allResTopIgnored,file="interimData/allResTopIgnored.RDat")
#alldf1TopIgnored=makeDF(bothMats,allResTopIgnored,motifFamilyTable)
#save(alldf1TopIgnored,file="interimData/alldf1TopIgnored.RDat")

load("interimData/allRes1Roadmap.RDat")
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemovedRoadmap.RDat")
alldf1=makeDF(bothMats,allRes1,motifFamilyTable)
save(alldf1,file="interimData/alldfRoadmap.RDat")
