#######################################################
## input: interimData/OverlappedCELs/.+
## output: interimData/ram_my_exprs.RDat

## explanation: aggregates CEL files and 
## uses rma-package: http://bmbolstad.com/misc/ComputeRMAFAQ/ComputeRMAFAQ.html
## RMA is the Robust Multichip Average. It consists of three steps: a background adjustment, quantile normalization (see the Bolstad et al reference) and finally summarization. Some references (currently published) for the RMA methodology are: Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15 

### only performs analysis for CEL-files that have a corresponding dhs-file
#######################################################
library(data.table)
library(oligo)
library(pd.huex.1.0.st.v2)
library(ff)
load("interimData/fileNameTables.RDat")

print("reading CELfiles")
celFiles <- list.celfiles('Data/CELs', full.names=TRUE)
celFilesWithoutPath=sub("Data/CELs/","",celFiles)

ind=is.element(celFilesWithoutPath,unlist(subset(celNames,hasDHS,CELfileName)))
celFilesTrunc=celFiles[ind]
### todo change the file 
affyExpressionFS <- read.celfiles(celFilesTrunc, pkgname="pd.huex.1.0.st.v2")


geneSummaries = rma(affyExpressionFS,  background=TRUE, normalize=TRUE,target="core") 
my_exprs=exprs(geneSummaries)
expressionAsRam=as.ram(my_exprs)
save(expressionAsRam, file='interimData/expressionAsRamDirect.RDat')
