library(mvtnorm)
library(MASS)
source("Code//fast_lmm.R")


fast_lmm_test1=function(){
###############################################################
###TEST1: test with no Covariance Structure ###################
###############################################################
set.seed(11)
print("TEST1: test with no Covariance Structure")
x=rnorm(70)
epsilon=rnorm(70)
y=x+epsilon*7
K=diag(70)
wer=fast_lmm(my_y=y,my_X=x,my_K=K)
wer2=summary(lm(y~x-1))
print("beta value of simple linear regression:")
print(wer2$coefficients[1])
print("beta value of mixed model regression:")
print(wer$beta)
print("pvalue value of simple linear regression:")
print(wer2$coefficients[4])
print("pvalue value of simple mixed model regression:")
print(1-pchisq(wer$chi_sq,1))
}

fast_lmm_test2=function(){
###############################################################
###TEST2: test with no Effect Size  ###########################
### but with Covariance structure  ############################
###############################################################
set.seed(11)
nrObs=250
nrOfTries=10000
toeplitzCoef=0.99
mixingCoef=0.5
toepSeq=rep(NA,nrObs)
for (i in c(0:(nrObs-1))){
    toepSeq[i+1]=toeplitzCoef^i
}
K=toeplitz(toepSeq)
cor_comp=mvrnorm(n=nrOfTries, rep(0,nrObs),K)
diagMat=matrix(0,dim(K)[1],dim(K)[2])
diag(diagMat)=1
unccor_comp=mvrnorm(n=nrOfTries, rep(0,nrObs),diagMat)
y=(1-mixingCoef)*cor_comp+mixingCoef*unccor_comp
x1=rnorm(nrObs,0,1)
all_pvals=rep(NA, nrOfTries)
for (i in c(1:nrOfTries)){
    print(i)
    wer=fast_lmm(y[i,],x1, K)
    all_pvals[i]=1-pchisq(wer$chi_sq,1)
}
print("number of pvalues below pval 0.05:")
print(mean(all_pvals<0.05))
return(all_pvals)
}

fast_lmm_test3=function(){
###############################################################
###TEST3: test with no Effect Size  ###########################
### but with Covariance structure diag not one  ###############
###############################################################
set.seed(11)
nrObs=250
nrOfTries=10000
toeplitzCoef=0.99
mixingCoef=0.5
toeplitzSeq=rep(NA,nrObs)
for (i in c(0:(nrObs-1))){
    toepSeq[i+1]=toeplitzCoef^i
}
K=toeplitz(toepSeq)
cor_comp=mvrnorm(n=nrOfTries, rep(0,nrObs),K)
diagMat=matrix(0,dim(K)[1],dim(K)[2])
diag(diagMat)=1
unccor_comp=mvrnorm(n=nrOfTries, rep(0,nrObs),diagMat)

y=5*((1-mixingCoef)*cor_comp+mixingCoef*unccor_comp)
x1=rnorm(nrObs,0,1)
all_pvals=rep(NA, nrOfTries)
for (i in c(1:nrOfTries)){
    print(i)
    wer=fast_lmm(y[i,],x1, K)
    all_pvals[i]=1-pchisq(wer$chi_sq,1)
}
print("number of pvalues below pval 0.05:")
print(mean(all_pvals<0.05))
return(all_pvals)
}

fast_lmm_test1()
a=fast_lmm_test2()
b=fast_lmm_test3()
