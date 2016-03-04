naive_regression <- function(my_y, my_X){

    len=length(my_X[1,])
    pVals=rep(NA,len)
    for (i in c(1:len)){
        print(i)
        wer=summary(lm(my_y~my_X[,i]))
        x=1-pt(wer$coef[2,3],107)
        pVals[i]=x
    }
    return(pVals)
    
}
