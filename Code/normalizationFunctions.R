
plotDifferencesInCor=function(mat1,mat2){
if(sum(dim(mat1)!=dim(mat2))!=0){
        print("warning:mats not equal size") 
}
d=rep(0,10)
for(i in c(1:length(mat1[,1]))){
d[i]=(cor(mat1[i,],mat2[i,]))}
hist(d,30)
print("median correlation")
print(median(d))
}

getFirstEig=function(mat,i){
        ### note: scaling has to happen somewhere else
        if(i<1){
            print("warning: cannot remove less than 1 eigen vects")
            return(NULL)    
        }
#   browser()
        dims=dim(mat)
##        matDeMeaned=t(scale(t(mat),center=T,scale=F))
        eig=eigen(t(mat)%*%(mat)/length(mat[,1]))
        eig$values
        return(eig$vectors[,c(1:i)])
}

getFirstEigCov=function(mat,i){
        ### note: scaling has to happen somewhere else
        if(i<1){
            print("warning: cannot remove less than 1 eigen vects")
            return(NULL)    
        }
        eig=eigen(cov(mat))
        #eig$values
        return(eig$vectors[,c(1:i)])
}

regressEigenVectOut=function(mat1,mat2){
## regresses significant columns of mat2 out of mat1 
        mat1RegressedOut=mat1
        for (i in c(1:length(mat1[1,]))){
        print(i)
           a=summary(lm(mat1[,i]~mat2))
              mySignificant=which(a$coefficients[,4]<0.01)-1  
                        browser()
         if (length(mySignificant)==0){
           mat1RegressedOut[,i]=mat1[,i]
           }else{
  
                a=lm(mat1[,i]~mat2[,mySignificant])
                        mat1RegressedOut[,i]=a$residuals

}}
return(mat1RegressedOut)
}

removeVects=function(mat,vects){
     ### note: scaling has to happen somewhere else 
   dims=dim(mat)
#    mat=t(scale(t(mat),center=T,scale=F))
    myResidMat=mat
        betas=mat%*%vects
        if(is.null(dim(vects))){
                return(myResidMat-kronecker(betas,t(vects)))
        }
        for(i in c(1:length(betas[1,]))){
              myResidMat=myResidMat-kronecker(betas[,i],t(vects[,i]))
        }
    return(myResidMat)
}

QQnormalizeMatCol=function(resultMat){
    dims=dim(resultMat)
    outMat=matrix(0,nrow=dims[1],ncol=dims[2])
    outMat=resultMat
    for (i in c(1:dims[2])){
    outMat[,i]=qqnorm(resultMat[,i],plot.it=F)$x
    }
    return(outMat)
}

QQnormalizeMatRow=function(resultMat){
    dims=dim(resultMat)
    outMat=matrix(0,nrow=dims[1],ncol=dims[2])
    outMat=resultMat
    for (i in c(1:dims[1])){
    outMat[i,]=qqnorm(resultMat[i,],plot.it=F)$x
    }
    return(outMat)
}


scalerRound=function(mat,n,scaleFlag=TRUE){
        for (i in c(1:n)){
                    mat=scale(t(scale(t(mat),scale=scaleFlag)),scale=scaleFlag)
                }
            return(mat)
    }


averageColQQnormalization=function(mat){
    mat2=apply(res,2,sort)
    aver=rowMeans(mat2)
    mat3=apply(mat,2,rank, ties.method="first")
    mat4=mat3
    for(i in c(1:length(mat3[1,]))){
        mat4[,i]=mat2[mat3[,i]]
    }
    return(mat4)
}
