#######################################################
##
## map GeneNames onto probesets:
## keep only 1to1 mappings
##
#######################################################
load("interimData/expressionAsRamDirect.RDat")
exprMatRaw=expressionAsRam
probesetGenenameMap=read.table("interimData/probesetGenenameMapDirect.tbl",sep="\t",header=F)
areDuplicated=(duplicated(probesetGenenameMap[,1]) | duplicated(probesetGenenameMap[,2]))
probesetGenenameMapUnique=probesetGenenameMap[!areDuplicated,]
myIds=as.numeric(rownames(exprMatRaw))
contained=is.element(myIds,probesetGenenameMapUnique[,1])
exprMatTrunc=exprMatRaw[contained,]
myIdsLeftover=as.numeric(rownames(exprMatTrunc))
myMatches=match(myIdsLeftover,probesetGenenameMapUnique[,1])

newNames=probesetGenenameMapUnique[myMatches,2]
rownames(exprMatTrunc)=newNames

expressionUnaveraged=exprMatTrunc
save(expressionUnaveraged, file="interimData/expressionUnaveragedDirect.RDat")

