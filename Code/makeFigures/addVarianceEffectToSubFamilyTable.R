library(svglite)
library(data.table)
set.seed(1)
load("interimData/overallAveragedExpressionDirect.RDat")
load("interimData/completeMotifFamilyTable.RData")
load("interimData/alldf1SubFam.RDat")

averagedExpression=averagedExpression[is.element(rownames(averagedExpression),motifFamilyTable[,geneSymbol]),]
allVar=log(apply(averagedExpression,1,var))
varNormed=qnorm(rank(allVar)/(length(allVar)+1))
##
unSubFams=alldf1SubFam[,subFamily]
len=length(unSubFams)
maxVar=rep(NA,len)
for(i in c(1:len)){
     genes=motifFamilyTable[is.element(subFamily,unSubFams[i]),geneSymbol]
     maxVar[i]=max(varNormed[is.element(names(varNormed),genes)])
}
maxVar=rep(NA,len)
for(i in c(1:len)){
     genes=motifFamilyTable[is.element(subFamily,unSubFams[i]),geneSymbol]
     maxVar[i]=max(varNormed[is.element(names(varNormed),genes)])
}

alldf1SubFam[,maxvar:=maxVar]
alldf1SubFamWithVariance=alldf1SubFam
save(alldf1SubFamWithVariance,file="interimData/alldf1SubFamWithVariance.RDat")

