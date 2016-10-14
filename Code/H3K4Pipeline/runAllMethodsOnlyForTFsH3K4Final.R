load("interimData/completeMotifFamilyTable.RData")
source("Code/normalizationFunctions.R")
library(multicore)
library(ggplot2)
library(data.table)


runRegression=function(bothMats,noCorrection=FALSE){
    my_K=cov(bothMats[[2]])
    source("Code/makeFigures/helperFunctionsFinal.R")
    bothMats[[2]]=bothMats[[2]][is.element(rownames(bothMats[[2]]),motifFamilyTable[,geneSymbol]),]
    indexSetList=makeSetListWingender(motifFamilyTable,bothMats,quote(subFamily))
    source("Code/fast_lmm_with_covariates.R")
    source("Code/fast_lmm.R")
    sets=sort(unique(unlist(indexSetList)))
    allYs=list()
    for(i in c(1:length(bothMats[[1]][,1]))){
        allYs[[i]]=cbind(bothMats[[1]][i,],bothMats[[3]][i,])
    }
    if(noCorrection==FALSE){
        fun=function(y){    
            res=fast_lmm_group(my_y=y[,1], my_X=t(bothMats[[2]][sets,]),my_K=my_K,covariates=y[,2])
            return(res)
        }
    }else{
        fun=function(y){    
            res=fast_lmm(my_y=y[,1], my_X=t(bothMats[[2]][sets,]),my_K=my_K)
            return(res)
        }
    }    
    allRes=mclapply(allYs,fun,mc.cores=25)
    return(allRes)
}

load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
bothMatsDHS=bothMats
load("interimData/bothMatDirectsH3K4DeterministicPCMotif1PCsRemoved.RDat")

matcher=match(colnames(bothMats[[1]]),colnames(bothMatsDHS[[1]]))
motifDHSpruned=bothMatsDHS[[1]][,matcher]
motifNamesInDHS=rownames(bothMatsDHS[[1]])
bothMats[[1]]=bothMats[[1]][is.element(rownames(bothMats[[1]]),motifNamesInDHS),]
matcherrow=match(rownames(bothMats[[1]]),rownames(bothMatsDHS[[1]]))
ww=motifDHSpruned[matcherrow,]
ww=t(scale(t(ww)))
topVects=getFirstEig(ww,1)
ww=removeVects(ww,topVects)

bothMats[[3]]=QQnormalizeMatRow(ww)

save(bothMats,file="interimData/bothMatDirectsH3K4out.RDat")
allResH3K4=runRegression(bothMats,noCorrection=FALSE)
save(allResH3K4,file="interimData/allResH3K4.RDat")

allResH3K4_uncor=runRegression(bothMats,noCorrection=TRUE)
save(allResH3K4_uncor,file="interimData/allResH3K4_uncor.RDat")

bothMatsRev=list()
bothMatsRev[[1]]=bothMats[[3]]
bothMatsRev[[3]]=bothMats[[1]]
bothMatsRev[[2]]=bothMats[[2]]

allResH3K4Rev=runRegression(bothMatsRev,noCorrection=FALSE)
save(allResH3K4Rev,file="interimData/allResH3K4Rev.RDat")
