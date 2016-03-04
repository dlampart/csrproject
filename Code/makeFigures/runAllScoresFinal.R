library(multicore)
library(data.table)
source("Code/fast_lmm.R")
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
allVals=list()
len=length(bothMats[[1]][,1])
motifNames=rownames(bothMats[[1]])
for (i in c(1:len)){
    allVals[[motifNames[i]]]=bothMats[[1]][i,]
}
my_X=t(bothMats[[2]])                                
fun=function(my_y){
    ee2=fast_lmm(my_y, my_X, my_K=NULL)
}
allRes=mclapply(allVals,fun,mc.cores=25)
save(allRes,file="interimData/completeRun.RDat")
