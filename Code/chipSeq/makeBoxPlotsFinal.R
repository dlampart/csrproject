library(svglite)
library(ggplot2)
library(data.table)
load("interimData/chipSeqEnrichment.RDat")
tbl=copy(plottingTable)
tbl1=plottingTable[variable=="enrichment shuffled",variable:="random TF"]
tbl2=plottingTable[variable=="enrichment",variable:="annotated TF"]
tbl=rbind(tbl1,tbl2)
load("interimData/chipSeqEnrichmentSubfamily.RDat")
tabl=plottingTable[variable=="enrichment",list(variable,value)]
tbl=rbind(tbl,tabl[,variable:="annotated Subfamily"])
svglite("PaperDocs/Images/chipSeqEnrichment.svg")
maxVal=max(tbl[,log2(value)])
minVal=min(tbl[,log2(value)])
p=ggplot(tbl,aes(x=variable,y=log2(value),fill=variable))+ylim(c(minVal-0.5,maxVal+0.5))+geom_boxplot(width=0.3)+ylab("enrichment (log2)")+xlab("")+theme_classic()+theme(axis.title=element_text(size=18),axis.text=element_text(size=18,angle=45,hjust=1),legend.position="none")
plot(p)
dev.off()

