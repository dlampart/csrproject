library(svglite)
library(ggplot2)
library(data.table)
source("Code/makeFigures/showPowerIncreasePerSubfamily_fun.R")
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")

load("interimData/alldf1SubFamH3K4_uncor.RDat")
alldf1_uncor=alldf1SubFam
load("interimData/alldf1SubFamH3K4.RDat")
alldf1_histone=alldf1SubFam
load("interimData/alldf1SubFamH3K4_rev.RDat")
alldf1_rev=alldf1SubFam
#######################################
## truncate to same nr of observations.
#######################################
overlappingSubfams=intersect(intersect(alldf1_uncor[,subFamily],alldf1_histone[,subFamily]),alldf1_rev[,subFamily])

alldf1_uncor=alldf1_uncor[is.element(subFamily,overlappingSubfams),]
alldf1_histone=alldf1_histone[is.element(subFamily,overlappingSubfams),]
alldf1_rev=alldf1_rev[is.element(subFamily,overlappingSubfams),]

labelList=c("histone","histone dhs corrected","dhs histone corrected")
tblList=list()
tblList[[1]]=alldf1_uncor
tblList[[2]]=alldf1_histone
tblList[[3]]=alldf1_rev
maxRank=getMaximalPossibleRank()
allResultsDF=makeAllResultsDF(tblList,labelList)
nrsPerRankTable=makeNrsPerRankTable(allResultsDF,labelList)
cumSumTable=makeCumSumTable(nrsPerRankTable,labelList)
b1=getMaximalPossibleRank()
a2=1
svglite("PaperDocs/Images/showPowerCumulativePerSubfamily_histone.svg")
p=ggplot(data=cumSumTable[cum<55,],aes(x=rank,y=probability,group=Method,colour=Method))+geom_line(size=3) + coord_cartesian(xlim=c(0,30))+geom_abline(slope=a2/b1)+geom_abline(slope=(10*a2/b1)) + scale_y_continuous(labels = scales::percent,expand=c(0,0))+ theme_classic() + ylim(0,0.45)+theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"),legend.title=element_text(size=18,face="bold"),legend.text=element_text(size=18))+theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 0.5))
plot(p)
dev.off()
ggsave("PaperDocs/Images/showPowerCumulativePerSubfamily_histone.png")
