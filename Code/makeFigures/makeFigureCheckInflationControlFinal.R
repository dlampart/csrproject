library(data.table)
library(ggplot2)
source("Code/naive_regression.R")
source("Code/fast_lmm.R")
load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
motifName="COE1"
ind=which(grepl(motifName,rownames(bothMats[[1]])))[1]

my_y=bothMats[[1]][ind,]
my_X=t(bothMats[[2]])
tt=naive_regression(my_y,my_X)

ee=fast_lmm(my_y, my_X, my_K=NULL)
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
ind=which(grepl(motifName,rownames(bothMats[[1]])))[1]

my_y=bothMats[[1]][ind,]
my_X=t(bothMats[[2]])
ee2=fast_lmm(my_y, my_X, my_K=NULL)

len=length(my_X[1,])
geneNames=colnames(my_X)

tabl=data.table(naive=tt,controlledlm0=1-pchisq(ee$chi_sq,1),controlledlm1=1-pchisq(ee2$chi_sq,1),theoretical=c(1:len)/(len+1),names=geneNames)
rankedTabl=tabl[,list(names,rank(naive),rank(controlledlm0),rank(controlledlm1))]
###############################
###############################
###############################

bonf=-log10(0.05/length(geneNames))
plotTabl=data.table(y=tabl[,-log10(sort(naive))],x=tabl[,-log10(sort(theoretical))],tabl[order(naive),names])
EBF1Table=plotTabl[V3=="EBF1",]
plotTabl=plotTabl[(x>0.7 | runif(len)>0.9) & (x>0.5 | runif(len)>0.3) & (x>0.3 | runif(len)>0.3) & V3!="EBF1",]
ggplot(data=plotTabl,aes(x=x,y=y))+geom_point(aes(y=y,x=x),shape=21,size=5,alpha=0.4,colour="purple4")+geom_point(data=EBF1Table,aes(y=y,x=x),size=8,shape=5)+xlim(0,4.3)+ylim(0,12)+geom_abline(slope=1,intercept=0)+geom_hline(yintercept=bonf,size=1,linetype=2,colour="red")+xlab("-log10(p-value) theoretical")+ylab("-log10(p-value) empirical")+theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16))+geom_text(x=0,y=bonf, label="Bonferroni cutoff", size=7, vjust=-0.5,hjust=0,fontface=6)
ggsave(paste("PaperDocs/Images/figControlCheckPanelA_",motifName,".svg",sep=""))
ggsave(paste("PaperDocs/Images/figControlCheckPanelA_",motifName,".png",sep=""))
###############################
###############################
###############################

bonf=-log10(0.05/length(geneNames))
plotTabl=data.table(y=tabl[,-log10(sort(controlledlm0))],x=tabl[,-log10(sort(theoretical))],tabl[order(controlledlm0),names])
EBF1Table=plotTabl[V3=="EBF1",]
plotTabl=plotTabl[(x>0.7 | runif(len)>0.9) & (x>0.5 | runif(len)>0.3) & (x>0.3 | runif(len)>0.3) & V3!="EBF1",]
ggplot(data=plotTabl,aes(x=x,y=y))+geom_point(aes(y=y,x=x),shape=21,size=5,alpha=0.4,colour="purple4")+geom_point(data=EBF1Table,aes(y=y,x=x),size=8,shape=5)+xlim(0,4.3)+ylim(0,12)+geom_abline(slope=1,intercept=0)+geom_hline(yintercept=bonf,size=1,linetype=2,colour="red")+xlab("-log10(p-value) theoretical")+ylab("-log10(p-value) empirical")+theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16))+geom_text(x=0,y=bonf, label="Bonferroni cutoff", size=7, vjust=-0.5,hjust=0,fontface=6)
ggsave(paste("PaperDocs/Images/figControlCheckPanelB_",motifName,".svg",sep=""))
ggsave(paste("PaperDocs/Images/figControlCheckPanelB_",motifName,".png",sep=""))
###############################
###############################
###############################

bonf=-log10(0.05/length(geneNames))
plotTabl=data.table(y=tabl[,-log10(sort(controlledlm1))],x=tabl[,-log10(sort(theoretical))],tabl[order(controlledlm1),names])
EBF1Table=plotTabl[V3=="EBF1",]
plotTabl=plotTabl[(x>0.7 | runif(len)>0.9) & (x>0.5 | runif(len)>0.3) & (x>0.3 | runif(len)>0.3) & V3!="EBF1",]
ggplot(data=plotTabl,aes(x=x,y=y))+geom_point(aes(y=y,x=x),shape=21,size=5,alpha=0.4,colour="purple4")+geom_point(data=EBF1Table,aes(y=y,x=x),size=8,shape=5)+xlim(0,4.3)+ylim(0,12)+geom_abline(slope=1,intercept=0)+geom_hline(yintercept=bonf,size=1,linetype=2,colour="red")+xlab("-log10(p-value) theoretical")+ylab("-log10(p-value) empirical")+theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16))+geom_text(x=0,y=bonf, label="Bonferroni cutoff", size=7, vjust=-0.5,hjust=0,fontface=6)
ggsave(paste("PaperDocs/Images/figControlCheckPanelC_",motifName,".svg",sep=""))
ggsave(paste("PaperDocs/Images/figControlCheckPanelC_",motifName,".png",sep=""))
