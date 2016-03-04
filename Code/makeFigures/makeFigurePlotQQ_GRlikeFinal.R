library(svglite)
library(data.table)
library(ggplot2)

source("Code/fast_lmm.R")
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
motifName="GCR"
ind=which(grepl(motifName,rownames(bothMats[[1]])))[1]

my_y=bothMats[[1]][ind,]
my_X=t(bothMats[[2]])
ee2=fast_lmm(my_y, my_X, my_K=NULL)


len=length(my_X[1,])
geneNames=colnames(my_X)


tabl=data.table(controlledlm1=1-pchisq(ee2$chi_sq,1),theoretical=c(1:len)/(len+1),names=geneNames)
numberOfPlotted=5
nameTable=tabl[rank(controlledlm1)<numberOfPlotted,]
nameTable=nameTable[order(controlledlm1),]

ww=c("NR3C1","NR3C2","AR","PGR")
ww2=nameTable[,names]

tabl2=tabl[order(controlledlm1),]
tabl2[,theoretical:=sort(theoretical)]
svglite(paste("PaperDocs/Images/qqPlotER-Like_",motifName,".svg",sep=""))
p=ggplot(data=tabl2, aes(x=-log10((theoretical)),y=-log10((controlledlm1))))+geom_point(size=3,colour="blue")+geom_point(data=tabl2[is.element(names,ww),], size=3,colour="red")+geom_abline(postion="identity")+geom_text(data=tabl2[is.element(names,ww2),],aes(label=names),hjust=1, vjust=-0.2,size=7)+xlab("-log10 p-values theoretical")+ylab("-log10 p-values empirical")+theme_classic()+theme(axis.title=element_text(size=18,face="bold"),axis.text=element_text(size=16))+ geom_text(data=tabl2[is.element(names,ww),],aes(label=names),hjust=1, vjust=-0.2,size=7)
plot(p)
dev.off()
