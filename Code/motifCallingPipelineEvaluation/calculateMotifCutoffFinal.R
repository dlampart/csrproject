library(multicore)
library(ggplot2)
source("Code/chipSeq/processNames.R")
load("interimData/fileNameTables.RDat")
wer3=NULL
nameTable=rbind(haibNames,sydhNames)
len=length(nameTable[,motifIndex])
motifFiles=paste("interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/",nameTable[,motifNames],".sorted.bed",sep="")
chipSeqFiles=paste("Data/chipSeq/",nameTable[,lab],"/",nameTable[,sub(".gz","",fileName)],sep="")
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
count=1

calculateOverlap=function(currentRow){
    print(currentRow$motifFile)
    currentTmpFileName=tempfile()
    cmdMakeCurrentTmpFile=paste("cut -f1,2,3,4,6,7 ",currentRow$motifFile," > ",currentTmpFileName,sep="")
    system(cmdMakeCurrentTmpFile)
    cmdstr=paste("bedmap --ec --max --echo ",currentRow$dhsFile," ",currentTmpFileName, " | perl -ne '!/^NAN/ && print' | perl -ne '/(.+)\\|(chr\\d+\\t\\d+\\t\\d+\\t\\.\\t)/ && print \"$2$1\\n\"' | bedmap --ec --echo --indicator - ",currentRow$chipSeqFile, "| perl -ple 's/\\|/\\t/'",sep="")
    cat(paste(cmdstr,"\n"))
    dd=fread(cmdstr)
    if(is.null(dd)){return(NULL)}
    out=unlist(lapply(quantile(dd[,V5],c(1:100)/100),function(x){dd[V5>x,mean(V6)]}))
     system(paste("rm", currentTmpFileName))
    return(out)
}


calculateOverlapOnlyOneCutoff=function(currentRow){
    print(currentRow$motifFile)
    currentTmpFileName=tempfile()
    cmdMakeCurrentTmpFile=paste("cut -f1,2,3,4,6,7 ",currentRow$motifFile," > ",currentTmpFileName,sep="")
    system(cmdMakeCurrentTmpFile)
    cmdstr=paste("bedmap --ec --max --echo ",currentRow$dhsFile," ",currentTmpFileName, " | perl -ne '!/^NAN/ && print' | perl -ne '/(.+)\\|(chr\\d+\\t\\d+\\t\\d+\\t\\.\\t)/ && print \"$2$1\\n\"' | bedmap --ec --echo --indicator - ",currentRow$chipSeqFile, "| perl -ple 's/\\|/\\t/'",sep="")
    cat(paste(cmdstr,"\n"))
    dd=fread(cmdstr)
    if(is.null(dd)){return(NULL)}
    out=unlist(lapply(quantile(dd[,V5],c(1:100)/100),function(x){dd[V5>x,mean(V6)]}))
    ss=which(out>cutoff)
    if(length(which(out>cutoff))==0){
        return(NULL)
    }
    out=quantile(dd[,V5],which(out>cutoff)[1]/100)
    system(paste("rm", currentTmpFileName))
    return(out)
}

motifChipSeqPairList=prepareMotifChipSeqPairList(motifFiles,shuffleMotifFiles=F)
for(cutoff in c(0.35,0.5,0.7)){
    resCutoff=mclapply(motifChipSeqPairList,calculateOverlapOnlyOneCutoff,mc.cores=20)
    names(res)=nameTable[,motifNames]
    names(resCutoff)=nameTable[,motifNames]
    ww=names(unlist(resCutoff))
    aa=sub("\\..+","",ww)
    aainv=setdiff(names(resCutoff),aa)
    dd=sub("_..\\..+","",ww)
    bb=sub(".+\\.","",ww)
    AA=data.table(motifName=aa,cutoff=unlist(resCutoff))
    unique(sub("\\..+","",ww))
    write.table(AA, file=paste("interimData/chipSeqVals_",cutoff,".txt",sep=""),sep="\t",row.names=FALSE,col.names=TRUE)
}
