library(multicore)
library(data.table)
cutoffs=c(0.3,0.5,0.7)
chipseqvals=fread(paste("interimData/chipSeqValsFinal.txt",sep=""))
chipseqvals=chipseqvals[!is.na(cutoff0.7),]
chipseqvals=chipseqvals[,list(median(cutoff0.3),median(cutoff0.5),median(cutoff0.7)),by=motifName]
setnames(chipseqvals,c("motifName",paste("cutoff",cutoffs,sep="")))

for(my_cutoff in cutoffs){
load("interimData/fileNameTables.RDat")
#currentVals=ss[,max(cutoff),by=motifName]
currentVals=chipseqvals[,list(motifName,eval(parse(text=paste("cutoff",my_cutoff,sep=""))))]
dd=c(1:dim(currentVals)[1])
dhsPath="interimData/DHSnormalized/"
dhsFiles=list.files(dhsPath)
dhsFullPath=paste(dhsPath,dhsFiles,sep="")
dhsCellNames=sub("interimData/DHSnormalized/Normalized.","",dhsFullPath)

motifPath="interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/"
motifFiles=list.files(motifPath)
motifFiles=motifFiles[grepl("\\.sorted\\.bed$",motifFiles)]
motifFullPath=paste(motifPath,motifFiles,sep="")
motifCellNames=sub(motifPath,"",sub(".sorted.bed","",motifFullPath))

calculateOverlap=function(i){
    currentMotifFullPath=motifFullPath[i]
    print("motif:")
    tmpFileName=paste(currentMotifFullPath,"_overlaps.tmp",sep="")
    tmpFileName2=paste(currentMotifFullPath,"_overlaps_2.tmp",sep="")
    print(tmpFileName)
    system(paste("touch ",tmpFileName,sep=""))
    system(paste("rm ",tmpFileName,sep=""))    
    rr=currentVals[i,eval(parse(text=paste("cutoff",my_cutoff,sep="")))]
    print(i)    
    commandStr2=paste("more ",currentMotifFullPath,"  | cut -f1,2,3,4,5,6 | perl -ne '/\\t(\\d+\\.\\d+)$/ && ($1>",rr,") && print' > ",tmpFileName2,sep="")
    print("wer0")
    print(commandStr2)
    system(commandStr2)    
    for (j in c(1:length(dhsFullPath))){
        commandStr=paste("ExternalCode/bin/bedops --ec -e -1 ",dhsFullPath[j]," ",tmpFileName2," | wc -l >> ",tmpFileName,sep="")
        print(commandStr)
        system(commandStr)        
    }
    system(paste("rm ",tmpFileName2))        
    overlaps=read.table(tmpFileName)
    system(paste("rm ",tmpFileName,sep=""))
    return(overlaps$V1)
}
##################################################################
#### put values into matrix.
##################################################################
matched=match(currentVals[,motifName],motifCellNames)
motifCellNames=motifCellNames[matched]
motifFullPath=motifFullPath[matched]
#wer=mclapply(dd,calculateOverlap,mc.cores=5)
wer=lapply(dd,calculateOverlap)

completeMat=matrix(NA, nrow=length(dhsFullPath),ncol=length(motifFullPath))
for (i in c(1:length(wer))){
    completeMat[,i]=wer[[i]]
}
##################################################################
#### give correct cellline index.
##################################################################
ind=match(dhsCellNames,dhsNames$DHSfileName)
dhsCellIds=dhsNames$cellId[ind]
rownames(completeMat)=dhsCellIds
##################################################################
#### give correct motif index.
##################################################################
colnames(completeMat)=motifCellNames
##################################################################
unnormedMotifActivity=t(completeMat)
save(unnormedMotifActivity,file=paste("interimData//unnormedMotifActivityDirectChipSeq_",my_cutoff,".RDat",sep=""))
}
