#######################################################
## Code/mapResultsToTFClusters.R
##
## input: interimData/resultMatList.RDat
## input: Data/UniProtGeneSymbolTable.tbl (downloaded from http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22, columns: Entry; Entry name; Gene names(primary))
## input: interimData/tfClassifications.txt
##
## output: interimData/motifFamilyTable.RData
## output: interimData/completeMotifFamilyTable.RData
 
## explanation: maps motif-ids to geneIds via  uniprotIds
## (gives almost complete coverage). additionally,
## gives for each factor, to which subfamily and family it belongs.
## completeMotifFamilytTable.RData contains information about
## all tfs given in wingender; not just the ones that we have
## a motif for. 
#######################################################
## make interimData/motifFamilyTable.RData
#######################################################

library(data.table)
load("interimData/unnormedMotifActivityDirect.RDat")
tfNames=sub("_..","",rownames(unnormedMotifActivity))

uniprotGeneIDTable=fread("Data/UniProtGeneSymbolTable.tbl",header=T)
setnames(uniprotGeneIDTable,c("uniprotId","UniProtKB-ID","geneSymbol"))
uniprotGeneIDTable[,`UniProtKB-ID`:=sub("_HUMAN","",`UniProtKB-ID`)]

tfGeneIdTable=unique(subset(uniprotGeneIDTable,is.element(`UniProtKB-ID`,tfNames)))


motifClassificationTable=fread("interimData/tfClassifications.txt",sep="\t",header=F)
setnames(motifClassificationTable,c("V1","V2"),c("tfId","uniprotId"))
## remove ids of isoforms
motifClassificationTable=subset(motifClassificationTable,!grepl("(\\.\\d+){5}", tfId, perl=T))
motifClassificationTable[,subFamily:=sub("\\.\\d+$","", tfId, perl=T)]
motifClassificationTable[,family:=sub("\\.\\d+$","",subFamily, perl=T)]
setkey(motifClassificationTable,uniprotId)

geneId=subset(tfGeneIdTable,is.element(`uniprotId`,motifClassificationTable[,uniprotId]),)
setkey(tfGeneIdTable,uniprotId)
geneId=motifClassificationTable[tfGeneIdTable,]
motifFamilyTable=motifClassificationTable[tfGeneIdTable,]
save(motifFamilyTable,file="interimData/motifFamilyTable.RData")
#######################################################
## make interimData/completeMotifFamilyTable.RData
#######################################################

library(data.table)
load("interimData/unnormedMotifActivityDirect.RDat")
tfNames=sub("_..","",rownames(unnormedMotifActivity))

uniprotGeneIDTable=fread("Data/UniProtGeneSymbolTable.tbl",header=T)
setnames(uniprotGeneIDTable,c("uniprotId","UniProtKB-ID","geneSymbol"))
uniprotGeneIDTable[,`UniProtKB-ID`:=sub("_HUMAN","",`UniProtKB-ID`)]

#tfGeneIdTable=unique(subset(uniprotGeneIDTable,is.element(`UniProtKB-ID`,tfNames)))
tfGeneIdTable=unique(uniprotGeneIDTable)

motifClassificationTable=fread("interimData/tfClassifications.txt",sep="\t",header=F)
setnames(motifClassificationTable,c("V1","V2"),c("tfId","uniprotId"))
## remove ids of isoforms
motifClassificationTable=subset(motifClassificationTable,!grepl("(\\.\\d+){5}", tfId, perl=T))
motifClassificationTable[,subFamily:=sub("\\.\\d+$","", tfId, perl=T)]
motifClassificationTable[,family:=sub("\\.\\d+$","",subFamily, perl=T)]
setkey(motifClassificationTable,uniprotId)

geneId=subset(tfGeneIdTable,is.element(`uniprotId`,motifClassificationTable[,uniprotId]),)
setkey(tfGeneIdTable,uniprotId)
#geneId=motifClassificationTable[tfGeneIdTable,]
motifFamilyTable=motifClassificationTable[geneId,]
save(motifFamilyTable,file="interimData/completeMotifFamilyTable.RData")
