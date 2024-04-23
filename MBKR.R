MBKR<-function(X,Y){
  #INPUT
  #X: an n x p matrix
  #Y: an n x 1 vector
  #OUTPUT
  #MBKR correlation
  #UTILITY = int {corr(I(X <= x), I(Y <= y))}^2 dF(x)dF(y)
  n=dim(X)[1]
  p=dim(X)[2]
  UTILITY=matrix(0,p,1)
  Aone=matrix(rep(1,n),n,1)
  oneA=t(rep(1,n))
  indY=Y%*%oneA<=Aone%*%t(Y)
  FY=colMeans(indY)
  ctrY=indY-Aone%*%FY
  for (kk in 1:p){
    indX=X[,kk]%*%oneA<=Aone%*%t(X[,kk])
    FX=colMeans(indX)
    ctrX=indX-Aone%*%FX
    options(warn=-1) 
    Rcor=cor(indX,indY)
    Rcor= Rcor[-which(is.na(Rcor)==TRUE)]
    UTILITY[kk,]=mean(Rcor^2)
  }
  return(UTILITY)  
}
###Example
n=200;
p=5;
mu=rep(0,p);
rho=0.5;   
Sigma=rho^(abs(matrix(rep(c(1:p),each=p),ncol=p,byrow=T)
               -matrix(rep(c(1:p),each=p),ncol=p)))
library(mvtnorm)
X=rmvnorm(n=n,mean=mu, sigma=Sigma) #n*p
beta=rep(0,p);beta[1:4]=c(5,5,5,-15*sqrt(rho))
beta=as.matrix(beta)
Y=X%*%beta+rnorm(n,0,1)
UTILITY=MBKR(X,Y)

