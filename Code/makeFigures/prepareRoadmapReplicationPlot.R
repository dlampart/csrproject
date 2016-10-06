library(svglite)
library(ggplot2)
library(data.table)
set.seed(1)
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")
load("interimData/alldfSubFamRoadmap.RDat")
alldf1SubFamRM=alldf1SubFam
load("interimData/alldf1SubFam.RDat")
alldfSubFam=alldf1SubFam
maxRank=getMaximalPossibleRank()

fun=function(cutoff, my_label, noinv=TRUE,df1,df2){
    my_sample=which(df1[,is.element(subFamily,df2[(rank<cutoff)==noinv,subFamily])])
    dfCSR=df1[my_sample,length(subFamily),by=rank]
    len=length(my_sample)
    dfCSR[,nrs:=V1]
    setkey(dfCSR,rank)
    dfCSR[,dist:=cumsum(nrs)/sum(nrs)]
    auglabel=paste(my_label, " (# ",len,")" ,sep="")
    dfCSR[,pioneers:=auglabel]
    plottingTable=rbind(data.table(rank=0,V1=0,nrs=0,pioneers=auglabel,dist=0),dfCSR)
    return(plottingTable)
}

g_legend=function(a.gplot){
    tmp_ggplot=gtable(ggplot_build(a.gplot))
    leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend = tmp$grobs[[leg]]
    return(legend)
}

plottingTable=fun(10,"rank < 10",noinv=TRUE,alldf1SubFam,alldf1SubFamRM)
plottingTable=rbind(plottingTable,fun(20,"rank < 20",noinv=TRUE,alldf1SubFam,alldf1SubFamRM))
plottingTable=rbind(plottingTable,fun(30,"rank < 30",noinv=TRUE,alldf1SubFam,alldf1SubFamRM))
plottingTable=rbind(plottingTable,fun(60,"rank < 60",noinv=TRUE,alldf1SubFam,alldf1SubFamRM))
plottingTable=rbind(plottingTable,fun(60,"rank >= 60",noinv=FALSE,alldf1SubFam,alldf1SubFamRM))

all_labels=unique(plottingTable[,pioneers])
#colors=c("000","000","000","000","000")
colors=paste("grey", floor(seq(from=1,to=70,length.out=5)),sep="")
p=ggplot()
p=p+geom_line(data=plottingTable,aes(x=rank,y=dist,group=pioneers,colour=pioneers),size=2,linetype=1) + scale_colour_manual(values=colors)
p=p + ylim(c(0,1)) + ylab("cumulative distribution") + xlab("subfamily rank score") + theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16),legend.title=element_blank(),legend.text=element_text(size=16),legend.position=c(0.75,0.25))+ geom_abline(slope=1/maxRank )+ theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
print(p)
#browser()
ggsave("PaperDocs/Images/encodeEnrichmentInRoadmapStrata.svg",dpi=300,width=5,height=4)
ggsave("PaperDocs/Images/encodeEnrichmentInRoadmapStrata.png",dpi=300,width=5,height=4)
print("pioneer enrichment plotted")

plottingTable=fun(10,"rank < 10",noinv=TRUE,alldf1SubFamRM,alldf1SubFam)
plottingTable=rbind(plottingTable,fun(10,"rank >= 10",noinv=FALSE,alldf1SubFamRM,alldf1SubFam))

all_labels=unique(plottingTable[,pioneers])
#colors=c("000","000","000","000","000")
colors=paste("grey", floor(seq(from=1,to=70,length.out=2)),sep="")
p=ggplot()
p=p+geom_line(data=plottingTable,aes(x=rank,y=dist,group=pioneers,colour=pioneers),size=2,linetype=1) + scale_colour_manual(values=colors)
p=p + ylim(c(0,1)) + ylab("cumulative distribution") + xlab("subfamily rank score") + theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16),legend.title=element_blank(),legend.text=element_text(size=16),legend.position=c(0.75,0.25))+ geom_abline(slope=1/maxRank )+ theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
print(p)
#browser()
ggsave("PaperDocs/Images/roadmapEnrichmentInEncodeStrata.svg",dpi=300,width=5,height=4)
ggsave("PaperDocs/Images/roadmapEnrichmentInEncodeStrata.png",dpi=300,width=5,height=4)
print("pioneer enrichment plotted")
