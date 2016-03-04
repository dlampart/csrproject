library(multicore)
library(ggplot2)
source("Code/chipSeq/processNames.R")
load("interimData/fileNameTables.RDat")

nameTable=haibNames
#nameTable=sydhNames

load("interimData/motifFamilyTable.RData")
###add flag
load("interimData/unnormedMotifActivityDirect.RDat")
motifNames=rownames(unnormedMotifActivity)
motifInds=match(motifFamilyTable[,`UniProtKB-ID`,],sub("_..$","",motifNames))
motifFamilyTable[,motifs:=motifNames[motifInds]]
###END: add flag 

fun=function(x){
    len=length(x)
    out=rep("",len)
    if(len>1){
        for(i in c(1:len)){
            out[i]=sample(setdiff(x,x[i]),1)
        }
    }
    return(out)
}

perms=motifFamilyTable[,list(geneSymbol,fun(motifs)),by=subFamily]
tfInds=match(nameTable[,antibodyTrunc],perms[,geneSymbol])
nameTable[,shuffledMapping:=perms[tfInds,V2]]

#############################
nameTable=nameTable[shuffledMapping!="",]
nameTable[,motifNames:=shuffledMapping]
###########################

len=length(nameTable[,motifIndex])
motifFiles=paste("interimData/motifInstances/HOCOMOCOHIGHCONF/",nameTable[,motifNames],"_HighConf.bed",sep="")
chipSeqFiles=paste("Data/chipSeq/",nameTable[1,lab],"/",nameTable[,sub(".gz","",fileName)],sep="")
dhsFiles=paste("interimData/DHSnormalized/",sub("^wg","Normalized.wg",nameTable[,DHSfileName]),sep="")








prepareMotifChipSeqPairList=function(motifFiles,shuffleMotifFiles=F){    
    if(shuffleMotifFiles){
        motifFiles=sample(motifFiles)
    }
    motifChipSeqPairList=list()
    for(i in c(1:len)){
        motifChipSeqPairList[[i]]=data.table(motifFile=motifFiles[i],chipSeqFile=chipSeqFiles[i],dhsFile=dhsFiles[i],currentIndex=i)
    }
    return(motifChipSeqPairList)
}

calculateOverlap=function(currentRow){
    print("motif:")
    tmpFileName=paste(currentRow[1,currentIndex],"__1_tmp.txt",sep="")
    system(paste("touch ",tmpFileName,sep=""))
    system(paste("rm ",tmpFileName,sep=""))
    tmpFileName2=paste(currentRow[1,currentIndex],"__2_tmp.txt",sep="")
    system(paste("touch ",tmpFileName2,sep=""))
    system(paste("rm ",tmpFileName2,sep=""))
    commandStr=paste("bedops --ec -e -1 ",currentRow[1,dhsFile]," ",currentRow[1,chipSeqFile]," > ",tmpFileName,sep="")    
    system(commandStr)
    lenIn=0
    if(file.info(tmpFileName)$size>0){
        lenIn=dim(fread(tmpFileName))[1]
    }
    commandStr=paste("bedops --ec -e -1 ",tmpFileName," ",currentRow[1,motifFile]," > ",tmpFileName2,sep="")
    system(commandStr)
    lenInMotif=0
    if(file.info(tmpFileName2)$size>0){
        lenInMotif=dim(fread(tmpFileName2))[1]
    }
    commandStr=paste("bedops --ec -n -1 ",currentRow[1,dhsFile]," ",currentRow[1,chipSeqFile]," > ",tmpFileName,sep="")    
    system(commandStr)
    lenOut=0
    if(file.info(tmpFileName)$size>0){
        lenOut=dim(fread(tmpFileName))[1]
    }    
    commandStr=paste("bedops --ec -e -1 ",tmpFileName," ",currentRow[1,motifFile]," > ",tmpFileName2,sep="")
    system(commandStr)
    lenOutMotif=0
    if(file.info(tmpFileName2)$size>0){
        lenOutMotif=dim(fread(tmpFileName2))[1]
    }
    commandStr=paste("bedops --ec -e -1 ",currentRow[1,dhsFile]," ",currentRow[1,motifFile]," > ",tmpFileName,sep="")    
    system(commandStr)
    lenMotifTot=0
    if(file.info(tmpFileName)$size>0){
        lenMotifTot=dim(fread(tmpFileName))[1]            
    }        
    system(paste("rm ",tmpFileName,sep=""))        
    system(paste("rm ",tmpFileName2,sep=""))
    vect=c(lenIn,lenInMotif,lenOut,lenOutMotif,lenMotifTot)
    return(vect)
}

motifChipSeqPairList=prepareMotifChipSeqPairList(motifFiles,shuffleMotifFiles=F)
res=mclapply(motifChipSeqPairList,calculateOverlap,mc.cores=20)


###calculateOverlap(motifChipSeqPairList2[[13]])

tbl=NULL
for(i in c(1:length(res))){
    tbl=rbind(tbl,data.table(lenIn=res[[i]][1],lenInMotif=res[[i]][2],lenOut=res[[i]][3],lenOutMotif=res[[i]][4],lenMotifTot=res[[i]][5]))
}

myenrichment=tbl[,lenInMotif/lenIn]/tbl[,lenOutMotif/lenOut]
nameTable[,enrichment:=myenrichment]

motifChipSeqPairList=prepareMotifChipSeqPairList(motifFiles,shuffleMotifFiles=T)
res=mclapply(motifChipSeqPairList,calculateOverlap,mc.cores=20)
tbl=NULL
for(i in c(1:length(res))){
    tbl=rbind(tbl,data.table(lenIn=res[[i]][1],lenInMotif=res[[i]][2],lenOut=res[[i]][3],lenOutMotif=res[[i]][4],lenMotifTot=res[[i]][5]))
}

myenrichment=tbl[,lenInMotif/lenIn]/tbl[,lenOutMotif/lenOut]
nameTable[,`enrichment shuffled`:=myenrichment]

plottingTable=melt(nameTable[!is.na(`enrichment shuffled`) & !is.na(enrichment),list(`enrichment shuffled`,enrichment)])
save(plottingTable,file="interimData/chipSeqEnrichmentSubfamily.RDat")

