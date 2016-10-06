source("Code/makeFigures/runAllMethodsOnlyForTFs_fun.R")
load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
allResMinus=runRegression(bothMats,noCorrection=TRUE)
save(allResMinus,file="interimData/allResMinus.RDat")

load("interimData/bothMatDirectsDeterministicPCMotif0PCsRemoved.RDat")
allRes0=runRegression(bothMats,noCorrection=FALSE)
save(allRes0,file="interimData/allRes0.RDat")

load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemoved.RDat")
allRes1=runRegression(bothMats,noCorrection=FALSE)
save(allRes1,file="interimData/allRes1.RDat")
