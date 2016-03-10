library(svglite)
library(data.table)
library(ggplot2)
load("interimData/alldf1SubFamWithVariance.RDat")
topRanks=alldf1SubFamWithVariance[,rank==1]
labels=rep("All other subfamilies",length(topRanks))
labels[topRanks]="Top scored subfamilies"
alldf1SubFamWithVariance[,grouped:=labels]
ggplot(data=alldf1SubFamWithVariance, aes(x=))
 
ggplot(data=alldf1SubFamWithVariance,aes(x=maxvar,group=grouped,fill=grouped)) +geom_density(,alpha=0.5) + theme_classic()+theme(legend.title=element_blank(),legend.text=element_text(size=16),legend.position=c(0.25,0.8)) + xlab("variance (normalized)") + guides(fill = guide_legend(override.aes = list(colour = NULL)))+theme(legend.key = element_rect(colour = "black"),axis.title = element_text(size=18))
ggsave("PaperDocs/Images/showVarianceEffect.svg")
ggsave("PaperDocs/Images/showVarianceEffect.png")
printingTable=alldf1SubFamWithVariance[,list(subFamily,rank, topGene)]
setnames(printingTable,c("subFamily","CSR rank (subfamily level)","top gene in subfamily"))
write.table(printingTable,file="PaperDocs/supplementaryTable1.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)




