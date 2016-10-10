preparePanelData=function(alltablInput,varNames){
    tot=NULL
    for(i in c(1:length(varNames))){
        for(j in c(1:length(varNames))){
            stratStr=varNames[i]
            carStr=varNames[j]
            alltabl=copy(alltablInput)
            alltabl=alltabl[!is.na(eval(parse(text=stratStr))),]
            alltabl=alltabl[!is.na(eval(parse(text=carStr))),]
            stratifyingSubfams=alltabl[eval(parse(text=stratStr))<10,subFamily]
            print(i)
            print(j)        
            alltabl[,isInGroup:=is.element(subFamily,stratifyingSubfams)]
            maxRank=getMaximalPossibleRank()
            dfInGroup=alltabl[isInGroup==TRUE,length(subFamily),by=eval(parse(text=carStr))]
            dfInGroup[,rank:=parse]
            dfNotInGroup=alltabl[isInGroup!=TRUE,length(subFamily),by=eval(parse(text=carStr))]
            dfNotInGroup[,rank:=parse]
            dfInGroup[,nrs:=V1]
            dfNotInGroup[,nrs:=V1]
            dfInGroup[,stratifying:="upper stratum"]
            dfNotInGroup[,stratifying:="lower stratum"]
            setkey(dfInGroup,rank)
            setkey(dfNotInGroup,rank)
            dfInGroup[,dist:=cumsum(nrs)/sum(nrs)]
            dfNotInGroup[,dist:=cumsum(nrs)/sum(nrs)]
            plottingTable=rbind(dfInGroup,dfNotInGroup)
            plottingTable[,parse:=NULL]
            plottingTable[,strat:=stratStr]
            plottingTable[,car:=carStr]
            plottingTable=rbind(data.table(V1=0,rank=0,nrs=0,dist=0,stratifying="upper stratum",strat=stratStr,car=carStr),plottingTable)
            plottingTable=rbind(data.table(V1=0,rank=0,nrs=0,dist=0,stratifying="lower stratum",strat=stratStr,car=carStr),plottingTable)
            tot=rbind(plottingTable,tot)
        }
    }
    return(tot)
}
