library(svglite)

library(ggplot2)
library(data.table)
source("Code/makeFigures/showPowerIncreasePerSubfamily_fun.R")
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")
load("interimData/alldfMinusSubFam.RDat")
load("interimData/alldf1SubFam.RDat")
load("interimData/alldf0SubFam.RDat")

tblList=list()
tblList[[1]]=alldfMinusSubFam
tblList[[2]]=alldf1SubFam
tblList[[3]]=alldf0SubFam
labelList=c("not corrected","corrected","pc1 corrected")
maxRank=getMaximalPossibleRank()
allResultsDF=makeAllResultsDF(tblList,labelList)
nrsPerRankTable=makeNrsPerRankTable(allResultsDF,labelList)
cumSumTable=makeCumSumTable(nrsPerRankTable,labelList)
b1=getMaximalPossibleRank()
a2=1
svglite("PaperDocs/Images/showPowerCumulativePerSubfamily.svg")
p=ggplot(data=cumSumTable[cum<55,],aes(x=rank,y=probability,group=Method,colour=Method))+geom_line(size=3) + coord_cartesian(xlim=c(0,30))+geom_abline(slope=a2/b1)+geom_abline(slope=(10*a2/b1)) + scale_y_continuous(labels = scales::percent,expand=c(0,0))+ theme_classic() + theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),legend.title=element_text(size=18,face="bold"),legend.text=element_text(size=18))+theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.5))
plot(p)
dev.off()
ggsave("PaperDocs/Images/showPowerCumulativePerSubfamily.png")
