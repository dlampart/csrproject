source("Code/makeFigures/prepareRankingTablesPerSubFamily_fun.R")

load("interimData/allRes0.RDat")
allRes=allRes0
load("interimData/alldf0.RDat")
alldf=alldf0[geneNr=="first",]
alldf0SubFam=prepareAllDfSubFam(alldf,motifFamilyTable,allRes)
save(alldf0SubFam,file="interimData/alldf0SubFam.RDat")

load("interimData/allRes1.RDat")
allRes=allRes1
load("interimData/alldf1.RDat")
alldf=alldf1[geneNr=="first",]
alldf1SubFam=prepareAllDfSubFam(alldf,motifFamilyTable,allRes)
save(alldf1SubFam,file="interimData/alldf1SubFam.RDat")

load("interimData/allResMinus.RDat")
load("interimData/allRes1.RDat")
allRes=addFakeEntriesToMinus(allResMinus,rownames(allRes1[[1]]))
load("interimData/alldfMinus.RDat")
alldf=alldfMinus[geneNr=="first",]
alldfMinusSubFam=prepareAllDfSubFam(alldf,motifFamilyTable,allRes)
save(alldfMinusSubFam,file="interimData/alldfMinusSubFam.RDat")
