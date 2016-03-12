library(svglite)
library(ggplot2)
library(data.table)
load("interimData/chipSeqEnrichment.RDat")
tbl=copy(plottingTable)
labels=rep("Other motifs",length(plottingTable[,value]))
labels[plottingTable[,variable=="enrichment"]]="Correct TF motifs"
plottingTable[,variable:=labels]
tbl=copy(plottingTable)
load("interimData/chipSeqEnrichmentSubfamily.RDat")
tabl=plottingTable[variable=="enrichment",list(variable,value)]
tabl[,variable:="TF subfamily motifs"]
tblPrint=rbind(tbl[variable=="Correct TF motifs",],tabl)
tblPrint=rbind(tblPrint,tbl[variable=="Other motifs",])

svglite("PaperDocs/Images/chipSeqEnrichment.svg")
maxVal=max(tblPrint[,log2(value)])
minVal=min(tblPrint[,log2(value)])
p=ggplot(tblPrint,aes(x=variable,y=(value),fill=variable))+scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)))+geom_boxplot(width=0.3)+ylab("enrichment of ChIP bound motifs")+xlab("")+theme_classic()+theme(axis.title=element_text(size=18),axis.text.x=element_text(size=18,angle=45,hjust=1),axis.text.y=element_text(size=18),legend.position="none")+annotation_logticks(sides="l")
plot(p)
dev.off()
ggsave("PaperDocs/Images/chipSeqEnrichment.png")

