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


load("interimData/allResMinus.RDat")
load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
alldfMinus=makeDF(bothMats,allResMinus,motifFamilyTable)
save(alldfMinus,file="interimData/alldfMinus.RDat")
load("interimData/allRes0.RDat")
load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
alldf0=makeDF(bothMats,allRes0,motifFamilyTable)
save(alldf0,file="interimData/alldf0.RDat")
load("interimData/allRes1.RDat")
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
alldf1=makeDF(bothMats,allRes1,motifFamilyTable)
save(alldf1,file="interimData/alldf1.RDat")




