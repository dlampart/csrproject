
#######################################################
##Code/ExpressionProcessing/prepareProbesetToGeneSymbolMappingDirect.R

## input: Data/mapping_71.txt (Downloaded from http://www.stemformatics.org/contents/download_mappings)
## input: Data/martEnsg2hgnc.txt  from biomart
## output:interimData/probesetGenenameMapDirect.tbl

## explanation: maps the core-probeset-ids  to gene symbols.
## by intersecting mapping refseq-exon annotation to AffyExonProbesetCore Annotation.
#######################################################
library(data.table)
wer=fread("Data/mapping_71.txt")
setnames(wer,"V2","Ensembl Gene ID")
wer1=fread("Data/martEnsg2hgnc.txt")
setkeyv(wer,"Ensembl Gene ID")
setkeyv(wer1,"Ensembl Gene ID")
wer2=merge(wer,wer1)
wer2=wer2[`HGNC symbol`!="",]
wer2[,`Ensembl Gene ID`:=NULL]
wer3=unique(wer2)
###remove doubly mapped probes
#wer3[is.element(`HGNC symbol`,aa),]
aa=wer3[duplicated(V1),V1]
wer4=wer3[!is.element(V1,aa),]

aa=wer4[duplicated(`HGNC symbol`),`HGNC symbol`]
wer5=wer4[!is.element(`HGNC symbol`,aa),]
write.table(wer5,file="interimData/probesetGenenameMapDirect.tbl",sep="\t",quote=F,row.names=F,col.names=F)
