rm(list=ls())
source('R/library.R')
p <- 2
n <- 10000
w <- 200

# 50/50 for the first half, 75/25 for the rest
Y = c(
  matrix(c(2,1), 1, n/2), 
  matrix(c(2,1,1,1), 1, n/2) 
)


outGen1 <- genAbrupt(n, muW = n/2, p=p, seed=1)
outGen2 <- genAbrupt(n, muW = n/4, p=p, seed=2)
X <- matrix(NA,n,p)
X[Y==1,] <- outGen1$X[Y==1]
X[Y==2,] <- outGen2$X[Y==2]

paramsLDA <- list()
etas <- matrix(0,n,4)
chats <- rep(1,n)
scores <- rep(1,n)
for (i in 1:n){
  if (i > 1){
    outClass <- classify.LDA(X[i,], paramsLDA)  
    chats[i] <- outClass$chat
    scores[i] <- outClass$score
  }
  paramsLDA <- update.LDA(X[i,], Y[i], paramsLDA) 
  if (length(paramsLDA[[1]])>0){
    etas[i,1] <- paramsLDA[[1]]$eta  
  }
  if (length(paramsLDA[[2]])>0){
    etas[i,2] <- paramsLDA[[2]]$eta  
  }
  etas[i,3] <- paramsLDA[[3]]$eta
  etas[i,4] <- paramsLDA[[4]]$eta
  
}


par(mfrow=c(2,1))
plot(which(Y==1),X[Y==1,1],col='red',ylim=c(min(X),max(X)))
points(which(Y==2),X[Y==2,1],col='blue')

plot(smooth(etas[,1]), type='l',ylim=c(0,0.3),col='red')
lines(smooth(etas[,2]), col='blue')
lines(smooth(etas[,3]), col='green')
lines(smooth(etas[,4]),col='black')
legend('topleft', c('Class 1','Class 2','Cov','Prior'), col=c('red','blue','green','black'), lty=1)

print(mean(chats==Y))
