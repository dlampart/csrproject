source("Code/makeFigures/runAllMethodsOnlyForTFs_fun.R")
load("interimData/bothMatDirectsDeterministicPCMotif1PCsRemovedRoadmap.RDat")

allRes1=runRegression(bothMats,noCorrection=FALSE)
save(allRes1,file="interimData/allRes1Roadmap.RDat")
