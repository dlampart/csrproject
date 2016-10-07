library(multicore)
##load("interimData/fileNameTables.RDat")
system("mkdir interimData/roadmapDHS_sorted/")
dhsPath="Data/roadmapDHS_filtered/"
dhsNewPath="interimData/roadmapDHS_sorted/"
dhsFiles=list.files(dhsPath)
dhsFullPath=paste(dhsPath,dhsFiles,sep="")
for(el in dhsFullPath){
    cmdstr=paste("sort-bed ",el, " > ", sub(dhsPath,dhsNewPath,el),sep="")
    print(cmdstr)
    system(cmdstr)
}
dhsPath=dhsNewPath
dhsFiles=list.files(dhsPath)
dhsFullPath=paste(dhsPath,dhsFiles,sep="")
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
    for (j in c(1:length(dhsFullPath))){
        print(j)
        commandStr=paste("bedops --ec -e -1 ",dhsFullPath[j]," ",currentMotifFullPath," | wc -l >> ",tmpFileName,sep="")
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
completeMat=matrix(NA, nrow=length(dhsFullPath),ncol=length(motifFullPath))
for (i in c(1:length(wer))){
    completeMat[,i]=wer[[i]]
}

##################################################################
#### give correct cellline index.
##################################################################

rownames(completeMat)=sub("-DNase.imputed.narrowPeak.bed","",dhsFiles)
##################################################################
#### give correct motif index.
##################################################################
colnames(completeMat)=motifCellNames
##################################################################
unnormedMotifActivity=t(completeMat)
save(unnormedMotifActivity,file="interimData/unnormedMotifActivityRoadmapDirect.RDat")


