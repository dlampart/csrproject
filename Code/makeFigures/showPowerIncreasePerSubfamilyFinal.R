library(svglite)

library(ggplot2)
library(data.table)
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/alldfMinusSubFam.RDat")
alldfMinus=alldfMinusSubFam
load("interimData/alldf1SubFam.RDat")
alldf1=alldf1SubFam
load("interimData/alldf0SubFam.RDat")
alldf0=alldf0SubFam

############################
############################
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
allVals=allResultsDF[,list(rank/maxRank,rank,method)]
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

numberOfTestedSubFams = length(motifFamilyTable[is.element(geneSymbol,rownames(allRes[[1]])),unique(subFamily)])
numberOfSubFamsWithMotifs = length(allVals[method=="corrected",rank])

b1=DD4[,maxRank]
a1=DD4[,max(cum)]

DD4[,Method:=factor(method,levels=c("not corrected", "corrected", "pc1 corrected"))]
a2=1
DD4[,probability:=cum/max(cum)]

svglite("PaperDocs/Images/showPowerCumulativePerSubfamily.svg")
p=ggplot(data=DD4[cum<55,],aes(x=rank,y=probability,group=Method,colour=Method))+geom_line(size=3) + coord_cartesian(xlim=c(0,30))+geom_abline(slope=a2/b1)+geom_abline(slope=(10*a2/b1)) + scale_y_continuous(labels = scales::percent,expand=c(0,0))+ theme_classic() + theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),legend.title=element_text(size=18,face="bold"),legend.text=element_text(size=18))
plot(p)
dev.off()
ggsave("PaperDocs/Images/showPowerCumulativePerSubfamily.png")
