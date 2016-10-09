library(svglite)
library(ggplot2)
library(data.table)
set.seed(1)
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")

load("interimData/alldf1SubFam_55.RDat")
a=alldf1SubFam
setkey(a,subFamily)
load("interimData/alldf1SubFam_5.RDat")
b=alldf1SubFam
setkey(b,subFamily)
#alltabl0=merge(a,b,suffixes=c('.5e-5','.1e-5'))
alltabl0=merge(a,b,suffixes=c('.5e5','.1e5'))
load("interimData/alldf1SubFam_6.RDat")
b=alldf1SubFam
setkey(b,subFamily)
alltabl=merge(alltabl0,b)
#alltabl[,`rank.1e-6`:=rank]
alltabl[,`rank.1e6`:=rank]
cutoffs=paste("rank.",c("5e5","1e5","1e6"),sep="")
cutoffs=paste("rank.",c("5e5","1e5","1e6"),sep="")

tot=NULL 
for(i in c(1:3)){
    for(j in c(1:3)){
        stratStr=cutoffs[i]
        carStr=cutoffs[j]        
        stratifyingSubfams=alltabl[eval(parse(text=stratStr))<10,subFamily]
        print(i)
        print(j)        
        alltabl[,isInGroup:=is.element(subFamily,stratifyingSubfams)]
        maxRank=getMaximalPossibleRank()
        dfInGroup=alltabl[isInGroup==TRUE,length(subFamily),by=eval(parse(text=carStr))]
        dfInGroup[,rank:=parse]
        dfNotInGroup=alltabl[isInGroup!=TRUE,length(subFamily),by=eval(parse(text=carStr))]
        dfNotInGroup[,rank:=parse]
        dfInGroup[,nrs:=V1]
        dfNotInGroup[,nrs:=V1]
        dfInGroup[,stratifying:="upper stratum"]
        dfNotInGroup[,stratifying:="lower stratum"]
        setkey(dfInGroup,rank)
        setkey(dfNotInGroup,rank)
        dfInGroup[,dist:=cumsum(nrs)/sum(nrs)]
        dfNotInGroup[,dist:=cumsum(nrs)/sum(nrs)]
        plottingTable=rbind(dfInGroup,dfNotInGroup)
        plottingTable[,parse:=NULL]
        plottingTable[,strat:=stratStr]
        plottingTable[,car:=carStr]
        plottingTable=rbind(data.table(V1=0,rank=0,nrs=0,dist=0,stratifying="upper stratum",strat=stratStr,car=carStr),plottingTable)
        plottingTable=rbind(data.table(V1=0,rank=0,nrs=0,dist=0,stratifying="lower stratum",strat=stratStr,car=carStr),plottingTable)
        tot=rbind(plottingTable,tot)
    }
}
tot[,strat:=sub("rank.","",sub("e","E-",strat))]
tot[,car:=sub("rank.","",sub("e","E-",car))]
pl=ggplot(tot[stratifying=="upper stratum",],aes(x=rank,y=dist,group=stratifying))+geom_line(size=2,linetype=1,colour="black") + geom_line(data=tot[stratifying=="lower stratum",],aes(x=rank,y=dist),colour="darkgrey",size=3)  + ylim(c(0,1)) + ylab("probability / stratifying") + xlab("rank") + theme_bw() + theme(axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),legend.title=element_blank(),legend.text=element_text(size=10),legend.position=c(0.85,0.7))+ geom_abline(slope=1/maxRank)+facet_grid(strat~car)
print(pl)
ggsave(paste("PaperDocs/Images/cutoffOverlapPlot_All.svg",sep=""),dpi=300,width=4,height=4)
ggsave(paste("PaperDocs/Images/cutoffOverlapPlot_All.png",sep=""),dpi=300,width=4,height=4)
