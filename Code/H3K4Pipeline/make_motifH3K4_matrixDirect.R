library(multicore)
load("interimData/histoneData/histoneNames.RDat")
histonePath="interimData/histoneData/rawbeds/"
histoneCellNames=histoneNames$V4
histoneFullPath=paste(histonePath,histoneCellNames,sep="")
motifPath="interimData/motifInstances/HOCOMOCOHIGHCONF/"
motifFiles=list.files(motifPath)
motifFullPath=paste(motifPath,motifFiles,sep="")
motifCellNames=sub(motifPath,"",sub("_HighConf.bed","",motifFullPath))

calculateOverlap=function(currentMotifFullPath){    
    print("motif:")
    tmpFileName=paste(currentMotifFullPath,"_overlaps.tmp",sep="")
    print(tmpFileName)
    system(paste("touch ",tmpFileName,sep=""))
    system(paste("rm ",tmpFileName,sep=""))
    for (j in c(1:length(histoneFullPath))){
        print(j)
        commandStr=paste("bedops --ec -e -1 ",histoneFullPath[j]," ",currentMotifFullPath," | wc -l >> ",tmpFileName,sep="")
        system(commandStr)
    }
    overlaps=read.table(tmpFileName)
    system(paste("rm ",tmpFileName,sep=""))
    return(overlaps$V1)
}

##################################################################
#### put values into matrix.
##################################################################
wer=mclapply(as.list(motifFullPath),calculateOverlap,mc.cores=40)
completeMat=matrix(NA, nrow=length(histoneFullPath),ncol=length(motifFullPath))
for (i in c(1:length(wer))){
    completeMat[,i]=wer[[i]]
}

##################################################################
#### give correct cellline index.
##################################################################
histoneNames$cellId


ind=match(histoneCellNames,histoneNames$V4)
histoneCellIds=histoneNames$cellId[ind]
rownames(completeMat)=histoneCellIds
##################################################################
#### give correct motif index.
##################################################################
colnames(completeMat)=motifCellNames
##################################################################
unnormedMotifActivity=t(completeMat)

save(unnormedMotifActivity,file="interimData/unnormedH3K4MotifActivityDirect.RDat")

