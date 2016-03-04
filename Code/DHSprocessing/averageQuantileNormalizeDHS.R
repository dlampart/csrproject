#######################################################
## Code/DHSprocessing/averageQuantileNormalizeDHS.R
##
## input: Data/DHSfdr0.01/$1.bed
## output:nterimData/DHSnormalized/Normalized.$1.bed
##
## explanation: cuts all DHS-files such that 'cutoff'-top
## peaks are kept. Then quantile Normalizes signal
## where a given quantile is the mean of the quantiles of
## all underlying empirical distributions
#######################################################
library(data.table)
cutoff=90000
DHSfiles=list.files("Data/DHSfdr0.01/", pattern=".bed")

i=1
meanSignal=rep(0,cutoff)
for (i in c(1:length(DHSfiles))){
    print(DHSfiles[i])
    loadStr=paste("Data/DHSfdr0.01/",DHSfiles[i],sep="")
    wer=fread(loadStr)
    setkey(wer,V7)
    len=length(wer[,V1])
    wer2=wer[c((len-cutoff+1):len),]
    meanSignal=meanSignal+wer2[,V7]
}
meanSignal=meanSignal/length(DHSfiles)

for (i in c(1:length(DHSfiles))){
    print(DHSfiles[i])
    loadStr=paste("Data/DHSfdr0.01/",DHSfiles[i],sep="")
    wer=fread(loadStr)
    setkey(wer,V7)
    len=length(wer[,V1])
    wer2=wer[c((len-cutoff+1):len),]
    wer2[,V5:=meanSignal]
    wer2[,c(1:5),with=F]
    write.table(wer2,"interimData/tmp", quote=F,sep="\t",row.names=F,col.names=F)
    outFileStr=paste("interimData/DHSnormalized/Normalized.",DHSfiles[i],sep="")
    systemStr=paste("sort-bed interimData/tmp > ",outFileStr,sep="")
    system(systemStr)

}

