library(svglite)
library(ggplot2)
library(data.table)

source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")
load("interimData/alldfMinus.RDat")
load("interimData/alldf1.RDat")
load("interimData/alldf0.RDat")

############################
############################
## prepare compareMethods figures
makeAllResultsDF=function(alldfMinus,alldf0,alldf1){    
    corrlevels=c("not corrected","corrected","pc1 corrected")
    allResultsDF=alldfMinus[,method:="not corrected"]
    allResultsDF=rbind(allResultsDF,alldf0[,method:="corrected"])
    allResultsDF=rbind(allResultsDF,alldf1[,method:="pc1 corrected"])
    allResultsDF[,Method:=factor(method,levels=corrlevels)]
    return(allResultsDF)

}
maxRank=getMaximalPossibleRank()
allResultsDF=makeAllResultsDF(alldfMinus,alldf0,alldf1)
allVals=allResultsDF[geneNr=="first",list(rank/maxRank,rank,method)]
allVals[,Method:=factor(method,levels=c("not corrected","corrected","pc1 corrected"))]

AA=allVals[,length(V1),by=list(rank,method)]
setkey(AA,rank)
BB=data.table(data.frame(rank=c(0,0,0),method=c("not corrected","corrected","pc1 corrected"),V1=c(0,0,0)))
CC=rbind(BB,AA)

DD=CC[method=="corrected"]
DD[,cum:=CC[method=="corrected",cumsum(V1)]]
DD2=CC[method=="not corrected"]
DD2[,cum:=CC[method=="not corrected",cumsum(V1)]]
DD3=CC[method=="pc1 corrected"]
DD3[,cum:=CC[method=="pc1 corrected",cumsum(V1)]]
DD4=rbind(DD,DD2,DD3)

b1=DD4[,maxRank]
a1=DD4[,max(cum)]
print(a1)

a2=1
DD4[,probability:=cum/max(cum)]
DD4[,Method:=factor(method,levels=c("not corrected", "corrected", "pc1 corrected"))]
svglite(file = "PaperDocs/Images/showPowerCumulative.svg")
p=ggplot(data=DD4[cum <= 150],aes(x=rank,y=probability,group=Method,colour=Method))+geom_line(size=3)+coord_cartesian(xlim=c(0,50))+geom_abline(slope=a2/b1) + scale_y_continuous(labels = scales::percent,expand=c(0,0)) + geom_abline(slope=(10*a2/b1))+theme_classic() + theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),legend.title=element_text(size=18,face="bold"),legend.text=element_text(size=18))+ expand_limits(x = 0, y = 0)
plot(p)
dev.off()


