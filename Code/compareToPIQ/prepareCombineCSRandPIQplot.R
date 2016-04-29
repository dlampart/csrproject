library(svglite)
library(ggplot2)
library(data.table)
source("Code/makeFigures/helperFunctionsFinal.R")
pioneerSubFamilies=getZaretPioneerSubFamilies()

piqResultsTable=fread("Data/piqResultsTable.txt")
load("interimData/completeMotifFamilyTable.RData")
setkey(piqResultsTable,geneSymbol)
setkey(motifFamilyTable,geneSymbol)
merged=merge(piqResultsTable,motifFamilyTable)
subfamilyPIQtable=merged[,max(pioneerScore),by=subFamily]
setnames(subfamilyPIQtable,"V1","subfamilyPIQscore")
load("interimData/alldf1SubFam.RDat")
setkey(alldf1SubFam,subFamily)
setkey(subfamilyPIQtable,subFamily)
resultsSubFam=merge(alldf1SubFam, subfamilyPIQtable)
EE=resultsSubFam[,list(rank(rank,ties.method="min"),rank(-subfamilyPIQscore,ties.method="min"))]
resultsSubFam[,scoreCombined:=apply(EE,1,max)]
resultsSubFam[,isPioneer:=is.element(subFamily,pioneerSubFamilies)]

valsPIQ=sort(resultsSubFam[,subfamilyPIQscore],decreasing=TRUE)
valsCSR=sort(resultsSubFam[,rank],decreasing=FALSE)
valsCombined=sort(resultsSubFam[,scoreCombined],decreasing=FALSE)
precBoth=rep(0,length(valsPIQ))
recBoth=rep(0,length(valsPIQ))
precCSR=rep(0,length(valsPIQ))
recCSR=rep(0,length(valsPIQ))
precPIQ=rep(0,length(valsPIQ))
recPIQ=rep(0,length(valsPIQ))
for(i in c(1:length(valsPIQ))){
    precBoth[i]=mean(resultsSubFam[scoreCombined <= valsCombined[i],isPioneer])
    recBoth[i]=sum(resultsSubFam[scoreCombined <= valsCombined[i],isPioneer])/sum(resultsSubFam[,isPioneer])
    precCSR[i]=mean(resultsSubFam[rank <= valsCSR[i],isPioneer])
    recCSR[i]=sum(resultsSubFam[rank <= valsCSR[i],isPioneer])/sum(resultsSubFam[,isPioneer])
    precPIQ[i]=mean(resultsSubFam[subfamilyPIQscore >= valsPIQ[i],isPioneer])
    recPIQ[i]=sum(resultsSubFam[subfamilyPIQscore >= valsPIQ[i],isPioneer])/sum(resultsSubFam[,isPioneer])
}
TT=data.table(unique(cbind(precBoth,recBoth,precCSR,recCSR,precPIQ,recPIQ)))
comb=TT[,list(precBoth,recBoth,"Combined")]
setnames(comb,c("Precision","Recall","Method"))
csr=TT[,list(precCSR,recCSR,"CSR rank")]
setnames(csr, c("Precision","Recall","Method"))
piq=TT[,list(precPIQ,recPIQ,"PIQ pioneer score")]
setnames(piq, c("Precision","Recall","Method"))
tbl=rbind(comb,csr)
tbl=rbind(tbl,piq)
ggplot(tbl,aes(x=Recall,y=Precision,colour=Method))+geom_line(size=2,linetype=1)+theme_bw()+theme(legend.position = c(0.75, 0.75),legend.background = element_rect(size=0.2, colour ="black",linetype="solid"))
ggsave("PaperDocs/Images/combineCSRandPIQ.pdf",dpi=300,width=4,height=4)

