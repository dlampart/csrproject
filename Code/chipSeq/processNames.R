library(data.table)
library(reshape2)
load("interimData/completeMotifFamilyTable.RData")
load("interimData/fileNameTables.RDat")
wer=fread("Data/chipSeq/wgEncodeSydhTfbsPrepared.txt",sep="\t")
wer2=wer[c(7:(length(wer[,V2])-1)),list(V1,V2,V3)]
casted=data.table(dcast(wer2,V1~V2,value.var="V3"))
filesToUse=casted[grepl("narrowPeak.gz",fileName),list(antibody,cell,view,fileName)]


filesToUse[,cellUpper:=gsub("-","",toupper(cell))]
dhsNames[,cellUpper:=DHSnameToUpper]
setkey(dhsNames,cellUpper)
setkey(filesToUse,cellUpper)
wer2=merge(dhsNames,filesToUse)
wer2[,antibodyTrunc:=toupper(sub("-","",sub("_.+$","",antibody)))]
load("interimData/bothMatDirectsProbPCMotif1PCsRemovedExpr0PCsNotRemoved.RDat")
tfNames=sub("_..$","",rownames(bothMats[[1]]))
tfGeneSymbols=motifFamilyTable[match(tfNames,`UniProtKB-ID`),geneSymbol]
motifTabl=data.table(motifNames=rownames(bothMats[[1]]),antibodyTrunc=tfGeneSymbols,motifIndex=c(1:length(tfNames)))
setkey(motifTabl,antibodyTrunc)
setkey(wer2,antibodyTrunc)
wer3=merge(wer2,motifTabl)



metaDataFile="Data/chipSeq/wgEncodeSydhTfbsPrepared.txt"
peakNames="narrowPeak.gz"




fun=function(metaDataFile,peakNames){
    library(data.table)
    library(reshape2)
    load("interimData/completeMotifFamilyTable.RData")
    load("interimData/fileNameTables.RDat")
    wer=fread(metaDataFile,sep="\t")
wer2=wer[c(7:(length(wer[,V2])-1)),list(V1,V2,V3)]
    casted=data.table(dcast(wer2,V1~V2,value.var="V3"))
    filesToUse=casted[grepl(peakNames,fileName),list(antibody,cell,view,fileName)]
    filesToUse[,cellUpper:=gsub("-","",toupper(cell))]
    dhsNames[,cellUpper:=DHSnameToUpper]
    setkey(dhsNames,cellUpper)
    setkey(filesToUse,cellUpper)
    wer2=merge(dhsNames,filesToUse)
    wer2[,antibodyTrunc:=toupper(sub("-","",sub("_.+$","",antibody)))]
    load("interimData/bothMatDirectsProbPCMotif1PCsRemovedExpr0PCsNotRemoved.RDat")
    tfNames=sub("_..$","",rownames(bothMats[[1]]))
    tfGeneSymbols=motifFamilyTable[match(tfNames,`UniProtKB-ID`),geneSymbol]
    motifTabl=data.table(motifNames=rownames(bothMats[[1]]),antibodyTrunc=tfGeneSymbols,motifIndex=c(1:length(tfNames)))
    setkey(motifTabl,antibodyTrunc)
    setkey(wer2,antibodyTrunc)
    wer3=merge(wer2,motifTabl)
    if(metaDataFile=="Data/chipSeq/wgEncodeSydhTfbsPrepared.txt"){
        wer3[,lab:="Sydh"]
    }
    if(metaDataFile=="Data/chipSeq/wgEncodeHaibTfbsPrepared.txt"){
        wer3[,lab:="Haib"]
    }
    return(wer3)
}

metaDataFile="Data/chipSeq/wgEncodeSydhTfbsPrepared.txt"
peakNames="narrowPeak.gz"
sydhNames=fun(metaDataFile,peakNames)


metaDataFile="Data/chipSeq/wgEncodeHaibTfbsPrepared.txt"
peakNames="broadPeak.gz"
haibNames=fun(metaDataFile,peakNames)
#####################
#####################
#####################    




