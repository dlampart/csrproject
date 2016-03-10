library(data.table)
set.seed(1)
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")

addSubfamilyColumnToDf=function(alldf,motifFamilyTable){
    alldf[,tfName:=sub("_..$","",motifName)]
    tfSubFamilies=motifFamilyTable[match(alldf[,tfName],`UniProtKB-ID`),subFamily]
    alldf[,subFamily:=tfSubFamilies]
    return(alldf)
}

prepareAllResSubFam=function(alldfWithSubFamily,allRes){
    allRes=allRes[alldfWithSubFamily[,index]]
    alldf=alldfWithSubFamily
    subFams=alldfWithSubFamily[,subFamily]
    unSubFams=unique(subFams)
    allResSubFam=list()
    len=length(unSubFams)
    someMotifFromSameSubfamily=rep("",len)
    for(i in c(1:len)){
        inds=which(subFams==unSubFams[i])
        res=allRes[[inds[1]]]
        someMotifFromSameSubfamily[i]=alldf[inds[1],tfName]
        if(length(inds)>1){
            for(j in c(2:length(inds))){
                print(j)
                indMat=res < allRes[[inds[j]]]
                res[indMat]=allRes[[inds[j]]][indMat]
            }
        }
        allResSubFam[[i]]=res        
    }
    names(allResSubFam)=someMotifFromSameSubfamily
    return(allResSubFam)
}

addFakeEntriesToMinus=function(allResMinus,rowNames){
    len=length(allResMinus[[1]][,1])
    for(i in c(1:length(allResMinus))){
        allResMinus[[i]]$betas=rep(0,len)
        allResMinus[[i]]$deltas=rep(0,len)
        allResMinus[[i]]$all_sigmasqe=rep(0,len)
        allResMinus[[i]]=as.data.frame(allResMinus[[i]])
        rownames(allResMinus[[i]])=rowNames
    }
    return(allResMinus)
}

getAllDfSubFam=function(allResSubFam,motifFamilyTable){
    alldfSubFam=NULL
    len=length(allResSubFam)
    for(i in c(1:len)){
        print(i)
        ss=getRank(allResSubFam[[i]], names(allResSubFam)[i],motifFamilyTable)
        alldfSubFam=rbind(alldfSubFam,ss)
    }
    alldfSubFam$subFamily=motifFamilyTable[match(names(allResSubFam),`UniProtKB-ID`),subFamily]
    return(alldfSubFam)
}

prepareAllDfSubFam=function(alldf,motifFamilyTable,allRes){
    alldf=addSubfamilyColumnToDf(alldf,motifFamilyTable)
    allResSubFam=prepareAllResSubFam(alldfWithSubFamily=alldf,allRes)
    return(getAllDfSubFam(allResSubFam,motifFamilyTable))
}

load("interimData/allRes0.RDat")
allRes=allRes0
load("interimData/alldf0.RDat")
alldf=alldf0[geneNr=="first",]
alldf0SubFam=prepareAllDfSubFam(alldf,motifFamilyTable,allRes)
save(alldf0SubFam,file="interimData/alldf0SubFam.RDat")

load("interimData/allRes1.RDat")
allRes=allRes1
load("interimData/alldf1.RDat")
alldf=alldf1[geneNr=="first",]
alldf1SubFam=prepareAllDfSubFam(alldf,motifFamilyTable,allRes)
save(alldf1SubFam,file="interimData/alldf1SubFam.RDat")

load("interimData/allResMinus.RDat")
load("interimData/allRes1.RDat")
allRes=addFakeEntriesToMinus(allResMinus,rownames(allRes1[[1]]))
load("interimData/alldfMinus.RDat")
alldf=alldfMinus[geneNr=="first",]
alldfMinusSubFam=prepareAllDfSubFam(alldf,motifFamilyTable,allRes)
save(alldfMinusSubFam,file="interimData/alldfMinusSubFam.RDat")
