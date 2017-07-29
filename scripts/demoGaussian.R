rm(list=ls())
source('R/library.R')
source('R/library_gaussian.R')

p <- 5

n <- 5000
w <- 1000

outGen <- genAbrupt(n=n, muW=w, p=p, seed=2)
X <- outGen$X

params <- list()
etas <- rep(NA,n)
mus <- matrix(NA,n,p)
for (i in 1:n){
  params <- streamGaussian(X[i,], params, alpha=0.001) 
  etas[i] <- params$eta
  mus[i,] <- params$mu
}
par(mfrow=c(2,1))
plot(etas, type='l', main='Learning Rates', xlab = 'Time', ylab='Value')
plot(X[,1], type='l', main='Quantity Monitored', xlab='Time') 
lines(mus[,1],col='red')
