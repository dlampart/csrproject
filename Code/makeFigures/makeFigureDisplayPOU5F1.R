library(data.table)
library(ggplot2)
source("Code/naive_regression.R")
source("Code/fast_lmm_group_reg.R")
source("Code/fast_lmm.R")
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
motifName="PO5F"
ind=which(grepl(motifName,rownames(bothMats[[1]])))[1]

my_y=bothMats[[1]][ind,]
my_X=t(bothMats[[2]])
ind=which(grepl(motifName,rownames(bothMats[[1]])))[1]

my_y=bothMats[[1]][ind,]
my_X=t(bothMats[[2]])
ee2=fast_lmm(my_y, my_X, my_K=NULL)


len=length(my_X[1,])
geneNames=colnames(my_X)

tabl=data.table(controlledlm1=1-pchisq(ee2$chi_sq,1),theoretical=c(1:len)/(len+1),names=geneNames)
rankedTabl=tabl[,list(rank(controlledlm1))]

numberOfPlotted=4
bonf=-log10(0.05/length(geneNames))
textTable=tabl[rank(controlledlm1)<numberOfPlotted,]
textTable=textTable[order(controlledlm1),]
ggplot(data=tabl, aes(x=-log10(sort(theoretical)),y=-log10(sort(controlledlm1))))+geom_point(size=4,fill="black",shape=21)+geom_point(size=3,fill="white",shape=21)+geom_abline(postion="identity")+geom_text(data=textTable,aes(x=-log10(c(1:(numberOfPlotted-1))/length(geneNames)),y=-log10(controlledlm1),label=names),hjust=1, vjust=-0.2,size=7)+geom_hline(yintercept=bonf,size=1,linetype=2)+xlab("-log10(p-value) theoretical")+ylab("-log10(p-value) empirical")+geom_text(x=0,y=bonf, label="Bonferroni cutoff", size=7, vjust=-0.5,hjust=0,fontface=6)+theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16))
ggsave(paste("PaperDocs/Images/figure_",motifName,".png",sep=""))
ggsave(paste("PaperDocs/Images/figure_",motifName,".svg",sep=""))
