library(svglite)
library(ggplot2)
library(data.table)
set.seed(1)
source("Code/makeFigures/helperFunctionsFinal.R")
source("Code/motifCallingPipelineEvaluation/showOverlap_fun.R")
load("interimData/alldf1SubFam_0.3_chipseq.RDat")
a=alldf1SubFam
setkey(a,subFamily)
load("interimData/alldf1SubFam_0.7_chipseq.RDat")
b=alldf1SubFam
setkey(b,subFamily)
#alltabl0=merge(a,b,suffixes=c('.5e-5','.1e-5'))
alltabl0=merge(a,b,all=TRUE,suffixes=c('.spec0.3','.spec0.7'))
load("interimData/alldf1SubFam_0.5_chipseq.RDat")
b=alldf1SubFam
setkey(b,subFamily)
alltabl=merge(alltabl0,b, all=TRUE)
alltabl[,`rank.spec0.5`:=rank]
load("interimData/alldf1SubFam.RDat")
b=alldf1SubFam
setkey(b,subFamily)
alltabl=merge(alltabl,b, all=TRUE,suffixes=c('','.cutoff1e5'))
#alltabl[,`rank.1e-6`:=rank]
#alltabl[,`rank.cutoff1e5`:=rank]
#cutoffs=paste("rank.",c("5e-5","1e-5","1e-6"),sep="")
cutoffs=paste("rank.",c('spec0.3','spec0.5','spec0.7','cutoff1e5'),sep="")
source("Code/chipSeq/processNames.R")
load("interimData/fileNameTables.RDat")
wer3=NULL
nameTable=rbind(haibNames,sydhNames)
allmotifs=unique(nameTable[,motifNames])
protnames=unique(sub("_..$","",allmotifs))
load("interimData/completeMotifFamilyTable.RData")
subfams=motifFamilyTable[is.element(`UniProtKB-ID`,protnames),subFamily]
alltabl=alltabl[is.element(subFamily,subfams),]
tot=preparePanelData(alltabl, cutoffs)

rename=function(x){
    out=sub("spec","ChIP ", sub("rank.","",sub("1e"," 1E-",x)))
    return(out)
}
tot[,strat2:=rename(strat)]
tot[,car2:=rename(car)]
maxRank=getMaximalPossibleRank()
pl=ggplot(tot[stratifying=="upper stratum",],aes(x=rank,y=dist,group=stratifying))+geom_line(size=2,linetype=1,colour="black") + geom_line(data=tot[stratifying=="lower stratum",],aes(x=rank,y=dist),colour="darkgrey",size=3)  + ylim(c(0,1)) + ylab("probability / stratifying") + xlab("rank") + theme_bw() + theme(axis.title=element_text(size=10,face="bold"),axis.text=element_text(size=7),legend.title=element_blank(),legend.text=element_text(size=8),legend.position=c(0.85,0.7))+ geom_abline(slope=1/maxRank)+facet_grid(strat2~car2)+theme(strip.text.x = element_text(size = 8))
print(pl)
ggsave(paste("PaperDocs/Images/chipseqOverlapPlot_All.svg",sep=""),dpi=300,width=4,height=4)
ggsave(paste("PaperDocs/Images/chipseqOverlapPlot_All.png",sep=""),dpi=300,width=4,height=4)

alltablR=alltabl[!is.na(rank.cutoff1e5),]
m=sum(alltablR[,rank.cutoff1e5<10])
n=sum(alltablR[,!rank.cutoff1e5<10])
k=sum(alltablR[,!is.na(rank.spec0.7)])
q=sum(alltablR[,!is.na(rank.spec0.7) & rank.cutoff1e5<10])
pval0.7=1-phyper(q-1,m,n,k)
k=sum(alltablR[,!is.na(rank.spec0.3)])
q=sum(alltablR[,!is.na(rank.spec0.3) & rank.cutoff1e5<10])
pval0.3=1-phyper(q-1,m,n,k)
k=sum(alltablR[,!is.na(rank.spec0.5)])
q=sum(alltablR[,!is.na(rank.spec0.5) & rank.cutoff1e5<10])
pval0.5=1-phyper(q-1,m,n,k)
print("pval0.7")
print(pval0.7)
print("pval0.5")
print(pval0.5)
print("pval0.3")
print(pval0.3)
