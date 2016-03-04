library(data.table)
#######################################################
## get all CEL-file names and derive celline-specifier

## this file is mapping names manually. Couldn't find the
## meta-data file.
#######################################################
celfiles=list.files("Data/CELs", pattern=".CEL")
celPruned=celfiles
celPruned=sub("GSM\\d+_X_Hs_","",celPruned)
celPruned=sub("GSM\\d+_hg19_wgEncodeUwAffyExonArray","",celPruned)
celPruned=sub("GSM\\d+_E_AH3_068_01_","",celPruned)
celPruned=sub("GSM\\d+_","",celPruned)
celPruned=sub("RawDataRep\\d.CEL","",celPruned)
celPruned=sub("_E_.+.CEL","",celPruned)
celPruned=sub("_B\\d.CEL","",celPruned)
celPruned=sub("FibroP_AG","AG",celPruned)
celPruned=gsub("_","",celPruned)
celPruned=gsub("-","",celPruned)

celNames=data.table(CELfileName=celfiles,CELnameRaw=celPruned,CELnameToUpper=toupper(celPruned))
## check if toupper()  loses elements (is not injective)
qq=toupper(unique(celNames[,CELnameRaw]))
qq[duplicated(qq)]
celNames[CELnameToUpper==qq[duplicated(qq)],]
## result:yes maps LNCaP and LNCAP: => OK
celNames[,cellId:=CELnameToUpper]
#######################################################
## get all DHS-file names
#######################################################

dhsfiles=list.files("Data/DHSfdr0.01",pattern="bb.bed")
dhsPruned=dhsfiles
dhsPruned=sub("wgEncode(Uw|Duke|UWDuke)Dnase","",dhsPruned,perl=T)
dhsPruned=sub(".fdr01peaks.hg19.bb.bed","",dhsPruned,perl=T)
dhsNames=data.table(DHSfileName=dhsfiles, DHSnameRaw=dhsPruned,DHSnameToUpper=toupper(dhsPruned))
## check if toupper()  loses elements (is not injective)
qq=toupper(unique(dhsNames[,DHSnameRaw]))
qq[duplicated(qq)]
dhsNames[DHSnameToUpper==qq[duplicated(qq)],]
## result:no:  => OK
dhsNames[,cellId:=DHSnameToUpper]
setkey(dhsNames,cellId)
#######################################################
## change cellId for selected experiments where automatic
## mapping doesn't work.
#######################################################
#(SKNS vs SKNSHRA) confirmed
#(PANC vs PANC1) confirmed
#(GM69 vs GM06990) confirmed
#(WERI vs WERIRb1) confirmed
#(HRCE vs HRCEpiC) confirmed
#(Ntera2 vs NT2D1) confirmed

#unsure:: left out
#(MonocytesCD14RO01746 vs CD14)?
#(iPS vs  iPSCWRU1 iPSNIHi11 iPSNIHi7)?

celNames$cellId[celNames$cellId=="SKNS"]="SKNSHRA"
celNames$cellId[celNames$cellId=="PANC"]="PANC1"
celNames$cellId[celNames$cellId=="GM69"]="GM06990"
celNames$cellId[celNames$cellId=="WERI"]="WERIRB1"
celNames$cellId[celNames$cellId=="HRCE"]="HRCEPIC"
celNames$cellId[celNames$cellId=="NTERA2"]="NT2D1"

setkey(dhsNames,cellId)
setkey(celNames,cellId)

#######################################################
## CEL: get origin of plate from name
#######################################################
fromUW=grepl("_X_",celNames[,CELfileName])
origin=rep("",length(fromUW))
origin[fromUW]="Uw"
origin[!fromUW]="Duke"    
celNames[,origin:=origin]
#######################################################
## DHS: get origin of plate from name
#######################################################
origin=sub("wgEncode","",dhsNames[,DHSfileName])
origin=sub("Dnase.+bed","",origin)
dhsNames[,origin:=origin]
#######################################################
## check if overlap is ok and save:
#######################################################
namesIntersected=unique(intersect(celNames[,cellId],dhsNames[,cellId]))
inBothForDHS=is.element(dhsNames[,cellId],namesIntersected)
dhsNames[,hasCEL:=inBothForDHS]

inBothForCEL=is.element(celNames[,cellId],namesIntersected)
celNames[,hasDHS:=inBothForCEL]

#######################################################
## check if overlap is ok and save:
#######################################################
print(unique(intersect(celNames[,cellId],dhsNames[,cellId])))
save(dhsNames,celNames,file="interimData/fileNameTables.RDat")

