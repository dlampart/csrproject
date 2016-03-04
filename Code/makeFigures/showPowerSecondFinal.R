library(svglite)
library(ggplot2)
library(data.table)

source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")
load("interimData/alldf1.RDat")
load("interimData/alldf1TopRemoved.RDat")

alldfTopRemoved=rbind(alldf1[geneNr=="second",],alldf1TopRemoved)
myMethod=alldfTopRemoved[,geneNr]
myMethod[myMethod=="first"]="top regressed out"
myMethod[myMethod=="second"]="top removed"
alldfTopRemoved[,Method:=myMethod]
maxRank=getMaximalPossibleRankForSecond()##
allVals=alldfTopRemoved[,list(rank/maxRank,rank,Method)]

allVals[,Method:=factor(Method,levels=c("top regressed out","top removed"))]
AA=allVals[,length(V1),by=list(rank,Method)]
setkey(AA,rank)
BB=data.table(data.frame(rank=c(0,0),Method=c("top regressed out","top removed"),V1=c(0,0)))
CC=rbind(BB,AA)

DD=CC[Method=="top regressed out"]
DD[,cum:=CC[Method=="top regressed out",cumsum(V1)]]
DD2=CC[Method=="top removed"]
DD2[,cum:=CC[Method=="top removed",cumsum(V1)]]
DD4=rbind(DD,DD2)

b1=DD4[,maxRank]
a1=DD4[,max(cum)]
print(a1)

a2=1
DD4[,probability:=cum/max(cum)]

#DD4[,Method:=factor(method,levels=c("not corrected", "corrected", "pc1 corrected"))]
svglite(file = "PaperDocs/Images/showPowerSecondCumulative.svg",height=10,width=10)
p=ggplot(data=DD4[cum <= 150],aes(x=rank,y=probability,group=Method,colour=Method))+geom_line(size=3)+coord_cartesian(xlim=c(0,50))+geom_abline(slope=a2/b1) + scale_y_continuous(labels = scales::percent,expand=c(0,0)) + geom_abline(slope=(10*a2/b1))+theme_classic() + theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),legend.title=element_text(size=18,face="bold"),legend.text=element_text(size=18),legend.position=c(0.8,0.55))+ expand_limits(x = 0, y = 0)+xlab("motif rank score")+ylab("cumulative distribution")
plot(p)
dev.off()


