library(data.table)
load("interimData/expressionUnaveragedDirect.RDat")
load("interimData/fileNameTables.RDat")

fileIds=colnames(expressionUnaveraged)
geneNames=rownames(expressionUnaveraged)
relCel=celNames[hasDHS==T,]

cellIds=relCel[,unique(cellId)]
av=matrix(0,ncol=length(cellIds),nrow=length(geneNames))
rownames(av)=geneNames
colnames(av)=cellIds

for(i in c(1:length(cellIds))){
    print(i)
    subInd=relCel[cellIds[i]==cellId,is.element(fileIds,CELfileName)]
    bb=expressionUnaveraged[,subInd]
    if(sum(subInd)>1){
        av[,i]=rowMeans(bb)
    }else{
        av[,i]=bb
    }
}


averagedExpression=av
save(averagedExpression, file="interimData/overallAveragedExpressionDirect.RDat")
