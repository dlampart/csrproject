load("interimData/completeMotifFamilyTable.RData")
library(multicore)
library(ggplot2)
library(data.table)

funDirectRegression=function(y,bothMats=bothMats,sets=sets){
    len=length(sets)
    print(y[1])
    allChisq=rep(0,len)
    for(i in c(1:len)){
        allChisq[i]=qchisq(1-summary(lm(y~bothMats[[2]][sets[i],]))$coef[2,4],1)
    }
    res=data.table(chi_sq=allChisq)
    rownames(res)=rownames(bothMats[[2]][sets,])
    return(res)
}

runRegression=function(bothMats,noCorrection=FALSE){
    my_K=cov(bothMats[[2]])
    source("Code/makeFigures/helperFunctionsFinal.R")
    bothMats[[2]]=bothMats[[2]][is.element(rownames(bothMats[[2]]),motifFamilyTable[,geneSymbol]),]
    indexSetList=makeSetListWingender(motifFamilyTable,bothMats,quote(subFamily))
    source("Code/fast_lmm.R")
    sets=sort(unique(unlist(indexSetList)))
    allYs=list()
    for(i in c(1:length(bothMats[[1]][,1]))){
        allYs[[i]]=bothMats[[1]][i,]
    }
    fun=function(y){
        res=fast_lmm(my_y=y, my_X=t(bothMats[[2]][sets,]),my_K=my_K)
        return(res)
    }
    if(noCorrection==FALSE){
        allRes=mclapply(allYs,fun,mc.cores=25)
#                allRes=lapply(allYs,fun)
    }else{
        funDirectRegressionWrap=function(y){
            funDirectRegression(y,bothMats=bothMats,sets=sets)
        }
        allRes=mclapply(allYs,funDirectRegressionWrap,mc.cores=25)
    }
    return(allRes)
}
