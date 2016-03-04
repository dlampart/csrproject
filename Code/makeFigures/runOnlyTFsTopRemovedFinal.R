load("interimData/completeMotifFamilyTable.RData")
library(multicore)
library(ggplot2)
library(data.table)


load("interimData/alldf1.RDat")


#################################
## add TopGene as covariate
#################################

runRegressionWithCovariate=function(bothMats,noCorrection=FALSE,priorDF){
    my_K=cov(bothMats[[2]])
    source("Code/makeFigures/helperFunctionsFinal.R")
    bothMats[[2]]=bothMats[[2]][is.element(rownames(bothMats[[2]]),motifFamilyTable[,geneSymbol]),]
    indexSetList=makeSetListWingender(motifFamilyTable,bothMats,quote(subFamily))
#################################
    ##getting covariate
#################################
    geneNames=rownames(bothMats[[2]])
    bothMats[[1]]=bothMats[[1]][priorDF[,index],]
    topGenePosition=match(priorDF[,topGene],geneNames)
    for(i in c(1:length(motifMat[,1]))){
        covariate=t(t(bothMats[[2]][is.element(rownames(bothMats[[2]]),priorDF[i,topGene])]))
        datList[[i]]=list(motifMat[i,],covariate,topGenePosition[i])
    }
#################################
    ##end covariate
#################################
    myIndexList=list()
    for(i in c(1:length(sets))){
        myIndexList[[i]]=i
    }
    fun=function(dat){
        res=fast_lmm_group(my_y=dat[[1]], my_X=t(bothMats[[2]]),indexList=myIndexList[setdiff(sets,dat[[3]])],my_K=my_K,covariates=dat[[2]])
        return(res)
    }
    source("Code/fast_lmm_group.R")
    allResWithCovariate0=mclapply(datList,fun,mc.cores=25)
    allResWithCovariate=list()
    for(i in c(1:length(allResWithCovariate0))){
        res=as.data.frame(allResWithCovariate0[[i]])
        rnames=rownames(bothMats[[2]][setdiff(sets,datList[[i]][[3]]),])
        rownames(res)=rnames    
        allResWithCovariate[[i]]=res    
    }    
    return(allResWithCovariate)
}


allRes1TopRemoved=runRegressionWithCovariate(bothMats=bothMats,noCorrection=FALSE,priorDF=alldf1[geneNr=="first",])
save(allRes1TopRemoved,file="interimData/allRes1TopRemoved.RDat")

load("interimData/alldf1.RDat")

alldf1Shuffled=copy(alldf1)
allTopGenes=alldf1Shuffled[,topGene]
alldf1Shuffled[,topGene:=sample(allTopGenes)]
allRes1TopRemovedShuffled=runRegressionWithCovariate(bothMats=bothMats,noCorrection=FALSE,priorDF=alldf1Shuffled[geneNr=="first",])
save(allRes1TopRemovedShuffled,file="interimData/allRes1TopRemovedShuffled.RDat")
