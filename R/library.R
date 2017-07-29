update.LDA <- function(
  x.new, c.new, paramsLDA, alpha=0.001
){
  # currently only supporting two classes
  stopifnot(is.element(c.new,c(1,2)))
  
  p <- length(x.new)
  # initialise
  if (length(paramsLDA)==0){
    paramsLDA <- list(
      cone=list(),
      ctwo=list(),
      poolcov=list(),
      prior=list()
    )
  }
  
  # update prior
  paramsLDA$prior <- streamMultinomial(c.new, paramsLDA$prior, alpha = alpha)
  
  # update class-conditionals, holding covariance fixed to the pool estimate
  if (length(paramsLDA$poolcov)==0){
    covFix <- diag(p)  
  } else {
    covFix <- paramsLDA$poolcov$S
  }
  paramsLDA[[c.new]] <- streamGaussian(x.new, paramsLDA[[c.new]], alpha=alpha, covFix=covFix)
  
  # update pooled covariance, holding mean fixed to the resp. class-cond
  if (length(paramsLDA[[c.new]])==0){
    muFix <- rep(0,p)
  } else {
    muFix <- paramsLDA[[c.new]]$mu
  }
  paramsLDA$poolcov <- streamGaussian(X[i,], paramsLDA$poolcov, alpha=alpha/10, muFix=muFix)
  
  return(paramsLDA)
}



classify.LDA <- function(
  x.new, paramsLDA
){
  p <- length(x.new)
  Js <- rep(NA,2)
  
  ldS <- paramsLDA[[3]]$ldS
  invS <- paramsLDA[[3]]$invS
  stopifnot(!is.null(ldS))
  stopifnot(!is.null(invS))
  mu1 <- paramsLDA[[1]]$mu
  mu2 <- paramsLDA[[2]]$mu
  if (is.null(mu1)){mu1 <- rep(0,p)}
  if (is.null(mu2)){mu2 <- rep(0,p)}
  pmu <- paramsLDA[[4]]$mu
  
  J1 <- 0.5*(ldS + t(x.new - mu1) %*% invS %*% (x.new - mu1))
  J2 <- 0.5*(ldS + t(x.new - mu2) %*% invS %*% (x.new - mu2))
  score <- J1-log(pmu[1])-(J2-log(pmu[2]))
  if (score > 0){
    chat <- 2  
  } else {
    chat <- 1
  }
  return(list(chat=chat,score=score))
}

streamMultinomial <- function(
  c.new=c.new, my=list(), 
  eta.fixed=NA, eta.min=0.0001, eta.max=0.3, alpha=0.001
){
  k = 2 # we only work with binomial for now
  stopifnot(is.element(c.new,1:k))
  
  
  
  ## INITIALISATIONS if t=0
  if (length(my)==0){
    ## primary initialisations
    mu <- rep(1,k)/k
    J <- 0
    # gradients
    mug <- rep(0,k)
    Jg <- 0 
    
    ## if adaptive forgetting, lambda=lambda.max, else lambda=lambda.fixed 
    if (is.na(eta.fixed)){
      eta <- eta.min
    } else {
      eta <- eta.fixed
    }
    
    my <- list(
      mu = mu, J = J,
      mug=mug, Jg=Jg,
      eta=eta
    )
    
  }
  J <- 0
  Jg <- -sum(my$mug * (as.numeric(c.new==1:k)/my$mu-1));  
  
  mu <- (1-my$eta) * my$mu + my$eta * as.numeric(c.new==1:k)
  mug <- (1-my$eta) * my$mug + (as.numeric(c.new==1:k)-my$mu)
  
  
  
  # update lambda
  if (is.na(eta.fixed)){
    eta = my$eta - alpha*Jg;
  } else {eta = my$eta}
  eta = min(eta,eta.max);
  eta = max(eta,eta.min);
  
  # return updated parameters
  
  my <- list(
    mu = mu, J = J,
    mug=mug, Jg=Jg,
    eta=eta
  )
  
  return(my)
}


KL <- function(mu0,mu1,S0,S1){
  
  
  ## check dimensions
  if (!is.matrix(S1)) {
    S1 <- matrix(S1)
  }
  if (nrow(S1) != ncol(S1)) {
    stop(message = "S1 not square.\n")
  }
  if (!is.matrix(S0)) {
    S0 <- matrix(S0)
  }
  if (nrow(S0) != ncol(S0)) {
    stop(message = "S0 not square.\n")
  }
  if (any(dim(S0) != dim(S1))){
    stop(message = "S1 and S0 have different dimensions\n")
  }
  if (!is.matrix(mu1)){
    mu1 <- as.matrix(mu1)
  }
  if (!is.matrix(mu0)){
    mu0 <- as.matrix(mu0)
  }
  if (nrow(mu1)!=nrow(S1)){
    stop(message = "mu1 does not have the right dimensions (must be px1)\n")
  }
  if (nrow(mu0)!=nrow(S0)){
    stop(message = "mu1 does not have the right dimensions (must be px1)\n")
  }
  
  d <- nrow(S1)
  
  KL <- sum(diag(solve(S1)%*%S0)) + t(mu1-mu0) %*% solve(S1) %*% (mu1-mu0) - d - log(det(S0)/det(S1))
  
  
  return(KL)
}

genAbrupt <- function(
  n=1000, 
  p=2,
  muW=n/4, 
  covW=muW, 
  muSpread=10,
  covSpread=2*p, 
  seed=1
){
  if (!is.na(seed)){
    set.seed(seed)
  }
  if (is.na(covW)){covW <- n+1}
  if (is.na(muW)){muW <- n+1}
  mutrue <- matrix(NA,n,p)
  X <- matrix(NA,n,p)
  vartrue <- list()
  KLtrue <- rep(NA,n)
  
  for (i in 1:n){
    if (i%%covW == 1){
      vartrue[[i]] <- rWishart(1, covSpread, (1/covSpread)*diag(p))[,,1]     
    } else {
      vartrue[[i]] <- vartrue[[i-1]]
    }
    if (i%%muW == 1){
      mutrue[i,] <- rnorm(p,0,muSpread)  
    } else {
      mutrue[i,] <- mutrue[i-1,]      
    }
    if (i > 1){
      KLtrue[i] <- KL(mutrue[i,],mutrue[i-1,],vartrue[[i]],vartrue[[i-1]])  
    }
    X[i,] <- chol(vartrue[[i]]) %*% rnorm(p,0,1) + mutrue[i,]
     
  }
  #cat('Avg KL at change:',)
  groundTruth <- list(
    mutrue=mutrue, vartrue=vartrue, KLtrue=KLtrue
  )
  
  return(list(X=X,groundTruth=groundTruth))
}
