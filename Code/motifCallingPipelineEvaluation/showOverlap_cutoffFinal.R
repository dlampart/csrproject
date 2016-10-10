library(svglite)
library(ggplot2)
library(data.table)
set.seed(1)
source("Code/makeFigures/helperFunctionsFinal.R")
source("Code/motifCallingPipelineEvaluation/showOverlap_fun.R")
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
tot=preparePanelData(alltabl, cutoffs)

tot[,strat:=sub("rank.","",sub("e","E-",strat))]
tot[,car:=sub("rank.","",sub("e","E-",car))]
maxRank=getMaximalPossibleRank()
pl=ggplot(tot[stratifying=="upper stratum",],aes(x=rank,y=dist,group=stratifying))+geom_line(size=2,linetype=1,colour="black") + geom_line(data=tot[stratifying=="lower stratum",],aes(x=rank,y=dist),colour="darkgrey",size=3)  + ylim(c(0,1)) + ylab("probability / stratifying") + xlab("rank") + theme_bw() + theme(axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=10),legend.title=element_blank(),legend.text=element_text(size=10),legend.position=c(0.85,0.7))+ geom_abline(slope=1/maxRank)+facet_grid(strat~car)
print(pl)
ggsave(paste("PaperDocs/Images/cutoffOverlapPlot_All.svg",sep=""),dpi=300,width=4,height=4)
ggsave(paste("PaperDocs/Images/cutoffOverlapPlot_All.png",sep=""),dpi=300,width=4,height=4)
