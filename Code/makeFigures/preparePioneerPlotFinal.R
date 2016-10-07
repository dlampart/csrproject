library(svglite)
library(ggplot2)
library(data.table)
set.seed(1)
source("Code/makeFigures/helperFunctionsFinal.R")
load("interimData/completeMotifFamilyTable.RData")
load("interimData/alldf1SubFamWithVariance.RDat")
alldfSubFam=alldf1SubFamWithVariance
###############################
preparePioneerMotifTable=function(motifFamilyTable){
    allFactors=NULL
    FOXAgenes=motifFamilyTable[grepl("^3.3.1.1$",subFamily),]
    POU5genes=motifFamilyTable[grepl("^3.1.10.5$",subFamily),]
    SOXBgenes=motifFamilyTable[grepl("^4.1.1.2$",subFamily),]
    KLFgenes=motifFamilyTable[grepl("^2.3.1.2$",subFamily),]     
    GATAgenes=motifFamilyTable[grepl("^2.2.1.1$",subFamily),]
    ASCLgenes=motifFamilyTable[grepl("^1.2.2.2$",subFamily),]
    PAXgenes=motifFamilyTable[grepl("^3.2.1.1$",subFamily),]
    SPIgenes=motifFamilyTable [grepl("^3.5.2.5$",subFamily),]
    P53genes=motifFamilyTable[grepl("^6.3.1.0$",subFamily),]
    CLOCKgenes=motifFamilyTable["1.2.5.2"==subFamily,]    
    allFactors=NULL
    allFactors=rbind(allFactors,FOXAgenes)
    allFactors=rbind(allFactors,POU5genes)
    allFactors=rbind(allFactors,SOXBgenes)
    allFactors=rbind(allFactors,KLFgenes)
    allFactors=rbind(allFactors,GATAgenes)
    allFactors=rbind(allFactors,ASCLgenes)
    allFactors=rbind(allFactors,PAXgenes)
    allFactors=rbind(allFactors,SPIgenes)
    allFactors=rbind(allFactors,P53genes)
    allFactors=rbind(allFactors,CLOCKgenes)
    return(allFactors)
}
###################################################
## END: prepare list of pioneer tfs from iwafuchi-Doi & Zaret 2014 
## output allFactors
###################################################
allFactors=preparePioneerMotifTable(motifFamilyTable)

alldfSubFam[,isPioneer:=is.element(subFamily,allFactors[,subFamily])]
maxRank=getMaximalPossibleRank()
dfPioneer=alldfSubFam[isPioneer==TRUE,length(subFamily),by=rank]
dfNotPioneer=alldfSubFam[isPioneer!=TRUE,length(subFamily),by=rank]
dfPioneer[,nrs:=V1]
dfNotPioneer[,nrs:=V1]
dfPioneer[,pioneers:="Pioneer"]
dfNotPioneer[,pioneers:="not Pioneer"]

setkey(dfPioneer,rank)
setkey(dfNotPioneer,rank)
dfPioneer[,dist:=cumsum(nrs)/sum(nrs)]
dfNotPioneer[,dist:=cumsum(nrs)/sum(nrs)]

plottingTable=rbind(dfPioneer,dfNotPioneer)
plottingTable=rbind(data.table(rank=0,V1=0,nrs=0,pioneers="Pioneer",dist=0),plottingTable)
plottingTable=rbind(data.table(rank=0,V1=0,nrs=0,pioneers="not Pioneer",dist=0),plottingTable)

ggplot(plottingTable[pioneers=="Pioneer",],aes(x=rank,y=dist,group=pioneers))+geom_line(size=2,linetype=1,colour="black") + geom_line(data=plottingTable[pioneers=="not Pioneer",],aes(x=rank,y=dist),colour="darkgrey",size=3) + geom_point(data=plottingTable[pioneers=="Pioneer" & rank>0,],aes(x=rank,y=dist),colour="black",size=8) + geom_point(data=plottingTable[pioneers=="Pioneer" & rank>0,],aes(x=rank,y=dist),colour="white",size=5) + ylim(c(0,1)) + ylab("cumulative distribution") + xlab("subfamily rank score") + theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16),legend.title=element_blank(),legend.text=element_text(size=20),legend.position=c(0.85,0.7))+ geom_abline(slope=1/maxRank) + theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

ggsave("PaperDocs/Images/pioneerEnrichment.svg",dpi=300,width=4,height=4)
ggsave("PaperDocs/Images/pioneerEnrichment.png",dpi=300,width=4,height=4)
print("pioneer enrichment plotted")


quants=alldfSubFam[isPioneer==TRUE,quantile(maxvar,c(0:2)/2)]
MM=alldfSubFam[isPioneer!=TRUE & maxvar>quants[1] & maxvar<= quants[2],]
MM2=alldfSubFam[isPioneer!=TRUE & maxvar>quants[2] & maxvar<= quants[3],]

print("calculate p-values controlled")
################
gg=rep(0,50000)
gg2=rep(0,50000)
gg3=rep(0,50000)
gg4=rep(0,50000)
for(i in c(1:50000)){
gg[i]=sum(MM[sample(c(1:length(MM[,rank])),4),rank])
gg2[i]=sum(MM2[sample(c(1:length(MM2[,rank])),4),rank])
gg3[i]=sum(MM[sample(c(1:length(MM[,rank])),4),rank<10])
gg4[i]=sum(MM2[sample(c(1:length(MM2[,rank])),4),rank<10])
}

tt=sum(alldfSubFam[isPioneer==TRUE,rank])
tt2=sum(alldfSubFam[isPioneer==TRUE,rank<10])
print("sum rank")
print(mean(tt>=gg+gg2))
print("enrichment")
print(mean(tt2<=gg3+gg4))

print("calculate p-values not controlled")
MM3=alldfSubFam[isPioneer!=TRUE,]
################
gg=rep(0,50000)
gg2=rep(0,50000)
gg3=rep(0,50000)
gg4=rep(0,50000)
for(i in c(1:50000)){
gg[i]=sum(MM3[sample(c(1:length(MM3[,rank])),8),rank])
gg2[i]=sum(MM3[sample(c(1:length(MM3[,rank])),8),rank<10])
}
tt=sum(alldfSubFam[isPioneer==TRUE,rank])
tt2=sum(alldfSubFam[isPioneer==TRUE,rank<10])
print("sum rank")
print(mean(tt>=gg))
print("enrichment")
print(mean(tt2<=gg2))

