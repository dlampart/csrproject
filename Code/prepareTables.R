library(data.table)
## input: interimData/alldf1SubFam.RDat
## input: interimData/subfamilyNames.txt
## input: interimData/familyNames.txt
load("interimData/alldf1SubFam.RDat")
load("interimData/completeMotifFamilyTable.RDat")
alldf1SubFam=as.data.table(alldf1SubFam)
subfamilyNames=fread("interimData/subfamilyNames.txt",header=FALSE)
setnames(subfamilyNames,c("subfamilyId","subfamilyNames"))
familyNames=fread("interimData/familyNames.txt",header=FALSE)
setnames(familyNames,c("familyId","familyNames"))
subfamilyNames[,familyId:=sub("\\.\\d+$","",subfamilyId)]
####

missingSubfamilies=setdiff(alldf1SubFam[,subFamily],subfamilyNames[,subfamilyId])
ll=sub(".+\\.(\\d+)$","\\1",missingSubfamilies)
ll2=sub("(.+)\\.\\d+$","\\1",missingSubfamilies)
dd=familyNames[match(ll2,familyId),]
dd[,subfamilyId:=paste(familyId,ll,sep=".")]
dd[,subfamilyNames:=familyNames]
dd[,familyId:=sub("\\.\\d+$","",subfamilyId)]
subfamilyNames=unique(rbind(subfamilyNames,dd[,list(subfamilyId,subfamilyNames,familyId)]))
###
setkey(subfamilyNames,familyId)
setkey(familyNames,familyId)
nameTable=merge(familyNames,subfamilyNames)
setkey(nameTable,subfamilyId)
setnames(alldf1SubFam,"subFamily","subfamilyId")
setkey(alldf1SubFam,"subfamilyId")
resultsTable=merge(alldf1SubFam,nameTable)

genenamesInPaper=c("FOXA1","SPI1","GATA6","SOX2","OTX2","CEBPD","KLF4","STAT5A","MYOG","EBF1","P63","TFAP2C","RELB","POU5F1")
#FOXA: Pioneer Transcription Factors Target Partial DNA Motifs on Nucleosomes to Initiate Reprogramming.
#FOXA: Opening of compacted chromatin by early developmental transcription factors HNF3 (FoxA) and GATA-4.
#FOXA:  Pioneer transcription factors in cell reprogramming. Genes and Development.
#SPI1: PU.1 and C/EBPalpha/beta convert fibroblasts into macrophage-like cells.
#SPI1:  Pioneer transcription factors in cell reprogramming. Genes and Development.
#GATA6: Opening of compacted chromatin by early developmental transcription factors HNF3 (FoxA) and GATA-4.
#GATA: Pioneer transcription factors in cell reprogramming. Genes and Development.
#SOX2: Pioneer Transcription Factors Target Partial DNA Motifs on Nucleosomes to Initiate Reprogramming.
#SOX2: Pioneer transcription factors in cell reprogramming. Genes and Development.
#TFAP2C: Pioneer transcription factors in cell reprogramming. Genes and Development.
#TFAP2C:    AP-2γ regulates oestrogen receptor-mediated long-range chromatin interaction and gene transcription
#CEBPD:C/EBPbeta induces chromatin opening at a cell-type-specific enhancer.
#STAT: Common Docking Domain in Progesterone Receptor-B links DUSP6 and CK2 signaling to proliferative transcriptional programs in breast cancer cells
#OTX2: Reorganization of Enhancer Patterns in Transition from Naive to Primed Pluripotency
#EBF1: Pioneering Activity of the C-Terminal Domain of EBF1 Shapes the Chromatin Landscape for B Cell Programming
#P63: TP53 engagement with the genome occurs in distinct local chromatin environments via pioneer factor activity
#RELB: c-Rel: A pioneer in directing regulatory T-cell lineage commitment?
#NR3C1:Dynamic exchange at regulatory elements during chromatin remodeling underlies assisted loading mechanism.
#MYOG: Pioneer transcription factors in cell reprogramming
#POU5F1: Pioneer Transcription Factors Target Partial DNA Motifs on Nucleosomes to Initiate Reprogramming
#CLOCK:Pioneer transcription factors in cell reprogramming
discussedResultsTable=resultsTable[is.element(topGene,genenamesInPaper),]
setkey(discussedResultsTable,rank)
setkey(resultsTable,rank)
tableMain=discussedResultsTable[,list(subfamilyNames, topGene,rank)]
setnames(tableMain,c("Subfamily name","Top gene in subfamily","CAR rank (subfamily level)"))
write.table(tableMain,"PaperDocs/tableMain.txt",sep="\t",quote=F,row.names=FALSE,col.names=TRUE)
tableSupporting=resultsTable[,list(subfamilyNames, topGene,rank,subfamilyId)]
setnames(tableSupporting,c("Subfamily name","Top gene in subfamily","CAR rank (subfamily level)","Subfamily id"))
write.table(tableSupporting,"PaperDocs/supportingTable.txt",sep="\t",quote=F,row.names=FALSE,col.names=TRUE)
