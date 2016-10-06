### make it compliant with encode unnormedExpression.RDat file format
source("Code/normalizationFunctions.R")
library(moments)
library(data.table)
mart=fread("Data/martEnsg2hgnc.txt")
mart=unique(mart[`HGNC symbol`!="",list(`Ensembl Gene ID`,`HGNC symbol`)])
#expressions=read.table("Data/57epigenomes.RPKM.pc",header=T)
expressions=read.table("Data/57epigenomes.N.pc",header=T)
expressions$E000=NULL
gene_idsExpr=expressions$gene_id
tabl=data.table(`Ensembl Gene ID`=gene_idsExpr,ind=c(1:length(gene_idsExpr)))
setkey(tabl,`Ensembl Gene ID`)
setkey(mart,`Ensembl Gene ID`)
mapping=mart[tabl]
setkey(mapping,ind)
expressions$gene_id=NULL
expressions=as.matrix(expressions)
expressions=expressions[mapping[,!is.na(`HGNC symbol`)],]
hgnc=mapping[!is.na(`HGNC symbol`),`HGNC symbol`]
rownames(expressions)=hgnc
expressions=expressions[!duplicated(hgnc),]

expressionsLog=log(expressions[rowMeans(expressions) > 50,]+1)
averagedExpression=averageColQQnormalization(expressionsLog)
save(averagedExpression,file="interimData/unnormedExpressionsRoadmap.RDat")
