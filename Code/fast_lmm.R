######## fast_lmm ######
fast_lmm <- function(my_y, my_X, my_K=NULL){
##########
#######implements fast-lmm algorithm for population-structure correction.
#####arguments my_y: n*1 vector of phenotype values
### my_X: m*n matrix of predictor variable
### my_K: m*m relationship matrix 
#the  algorithm tries to estimate the scalar b where the my_y=my_x*b+epsilon
# and epsilon has distribution N(0, sigma_e^2*I+sigma_g^22*K)
if(is.null(my_K)){
    print("no relationship matrix given, using covariance of t(my_X)")
    my_K=cov(t(my_X))
}
    
##decompose my_K:
my_eig=eigen(my_K)
K_vectors=my_eig$vectors
K_values=my_eig$values
 


y_prime=t(K_vectors)%*%my_y
data=data.frame(y=y_prime, S=K_values)
likeli_h0=function(data,delta){	
	sigma_sq=(1/length(data$y))*sum(((data$y)^2)/(data$S+delta))
	return(length(data$y)*log(sigma_sq)+sum(log(data$S+delta))+sum(((data$y)^2)/(sigma_sq*(data$S+delta))))
}
lower_val=0.02
upper_val=10
optim_lm0=optim(par=0.5, fn=likeli_h0, data=data, method="Brent",lower=lower_val, upper=upper_val)
delta_lm0=optim_lm0$par
loglik_lm0=optim_lm0$value
if (delta_lm0==upper_val || delta_lm0==lower_val){
   print("warning delta_lm at edge of defined parameter_space")
   }
print(delta_lm0)
sigma_sq0=sigma_sq=(1/length(data$y))*sum(((data$y)^2)/(data$S+delta_lm0))
X_prime=t(K_vectors)%*%my_X

likeli_h=function(data,delta){	
	beta=(1/sum((data$x^2)/(data$S+delta)))*sum(data$x*data$y*(1/(data$S+delta)))
	sigma_sq=(1/length(data$y))*sum(((data$y-data$x*beta)^2)/(data$S+delta))
	loglik=(length(data$y)*log(sigma_sq)+sum(log(data$S+delta))+sum(((data$y-data$x*beta)^2)/(sigma_sq*(data$S+delta))))	
	return(loglik)	
}
my_n=dim(X_prime)[2]
all_betas=rep(0, my_n)
all_chisq=rep(0,my_n)
all_deltas=rep(0,my_n)
all_sigmasqe=rep(0,my_n)

for (i in c(1:my_n)){
if((i %% 100)==0){
    print(i)
}
data=data.frame(y=y_prime, x=X_prime[,i], S=K_values)
optim_lm=optim(par=0.5, fn=likeli_h, data=data, method="Brent",lower=lower_val, upper=upper_val)
delta_lm=optim_lm$par
loglik_lm=optim_lm$value
chisq=loglik_lm0-loglik_lm
pchisq(chisq,1)
#if (delta_lm==upper_val || delta_lm==lower_val){
#   print("warning delta_lm at edge of defined parameter_space")
#   }
beta_lm=(1/sum((data$x^2)/(data$S+delta_lm)))*sum(data$x*data$y*(1/(data$S+delta_lm)))
sigmasq_g_lm=(1/length(data$y))*sum(((data$y-data$x*beta_lm)^2)/(data$S+delta_lm))
sigmasq_e_lm=delta_lm*sigmasq_g_lm
all_betas[i]=beta_lm
all_chisq[i]=chisq
all_deltas[i]=delta_lm
all_sigmasqe[i]=sigmasq_e_lm
}

fast_lmm=data.frame(betas=all_betas,chi_sq=all_chisq,deltas=all_deltas,all_sigmasqe=all_sigmasqe)
rownames(fast_lmm)=colnames(my_X)

return(fast_lmm)
}
