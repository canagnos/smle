rm(list=ls())
source('R/library.R')
p <- 5
n <- 1000
w <- 200

outGen <- genAbrupt(n=n, muW=w, p=p, seed=2)
X <- outGen$X

params <- list()
etas <- rep(NA,n)
mus <- matrix(NA,n,p)
muFix <- rep(10,p)
covFix <- diag(p)
for (i in 1:n){
  params <- streamGaussian(X[i,], params, alpha=0.001, covFix=covFix) 
  etas[i] <- params$eta
  mus[i,] <- params$mu
}
par(mfrow=c(2,1))
plot(etas)
plot(X[,1]); lines(mus[,1],col='red')
