streamGaussian <- function(
  x.new=x.new, my=list(), 
  eta.fixed=NA, eta.min=0.001, eta.max=0.3, 
  alpha=0.001,
  lambda.fixed=NA, lambda.min=0.7, lambda.max=0.999,
  efficient=FALSE, scaleit=FALSE, 
  covFix=NULL, muFix=NULL # helper options for LDA
){
  
  ## x.new must be a vector
  ## params must be a list
  
  p = length(x.new); 
  
  ## INITIALISATIONS if t=0
  if (length(my)==0){
    
    
    
    ## primary initialisations
    v.0 <- 10
    nn = 1
    mu <- rep(0,p) 
    if (!is.null(muFix)){
      mu <- muFix
    } 
    S <- diag(v.0,p,p)
    if (!is.null(covFix)){
      S <- covFix
    }
    invS <- diag(1/v.0,p,p)
    ldS <- p*log(v.0) 
    J <- 0
    
    # gradients
    mug <- rep(0,p)
    Sg <- matrix(0,p,p)
    invSg <- matrix(0,p,p) 
    ldSg <- 0 
    Jg <- 0 
    
    ## if adaptive forgetting, lambda=lambda.max, else lambda=lambda.fixed 
    
    if (is.na(eta.fixed)){
      eta <- eta.min
      lambda <- lambda.max
    } else {
      eta <- eta.fixed
      lambda <- lambda.fixed
    }
    
    
    my <- list(
      mu = mu, S = S, invS = invS, ldS = ldS, J = J,
      mug=mug, Sg=Sg, invSg=invSg, ldSg=ldSg, Jg=Jg,
      eta=eta, lambda = lambda
    )
    
  }
  
  
  ## COMPUTE GRADIENT 
  J = my$ldS + t(x.new - my$mu)%*%my$invS%*%(x.new - my$mu)
  Jg = my$ldSg -2*t(x.new-my$mu)%*%(my$invS)%*%(my$mug) + t(x.new-my$mu)%*%my$invSg%*%(x.new-my$mu)
  if (scaleit) {Jg = (1/J)*Jg}
  
  ## PARAMETER UPDATES
  if (!is.null(muFix)){
    mug = rep(0,p)
    mu = muFix
  } else {
    mu = my$mu + my$eta*(x.new - my$mu)
    mug = (1-my$eta)*my$mug +(x.new-my$mu);    
  }
  
  if (!is.null(covFix)){
    S = covFix
    Sg = 0*diag(p)
  } else {
    S = (1-my$eta)*my$S + my$eta*((x.new-mu)%*%t(x.new-mu))
    Sg = (1-my$eta)*my$Sg -my$S + (x.new-mu)%*%t(x.new-mu) + 
      my$eta*( -mug%*%t(x.new-mu) - (x.new-mu)%*%t(mug))
    
    #Pi = my$Pi + my$eta*(x.new%*%t(x.new) - my$Pi) 
    #Pig = (1-my$eta)*my$Pig +(x.new%*%t(x.new)-my$Pi);
    
    #S = Pi - mu%*%t(mu)
    #Sg = Pig - mug%*%t(mu) - mu%*%t(mug)
    
  }
  
  
  
  
  
  if (efficient){
    ### NOT CHECKED YET !!!!!!
    # auxillary quantities gamma, h and their gradients
    g = (n-1)^2/n + t(x.new-mu)*my$invS*(x.new-mu)
    gg = ng*(n-1)*(n+1)/n^2+t(x.new-mu)*my$invSg*(x.new-mu)-2*t(x.new-mu)*my$invS*mug;
    
    h = my$invS*(x.new-mu)*t(x.new-mu)*my$invS;
    hg = my$invSg*(x.new-mu)*t(x.new-mu)*my$invS+my$invS*(x.new-mu)*t(x.new-mu)*my$invSg;
    hg = hg - my$invS*mug*t(x.new-mu)*my$invS- my$invS*(x.new-mu)*t(mug)*my$invS;
    
    # inv(sigma) and log(det(sigma)), updated efficiently
    invSg = -(ng/(n-1)^2)*(my$invS-(1/g)*h)+(n/(n-1))*(my$invSg+(gg/g^2)*h-(1/g)*hg);
    invS = (n/(n-1))*(my$invS -(1/g)*h);
    
    ldS = (p-2)*log(n-1)+ (1-p)*log(n)+log(g)+my$ldS;
    ldSg = (p-2)*ng/(n-1) + (1-p)*ng/n + gg/g+my$ldSg;
  } else {
    # inv(sigma) and log(det(sigma)), updated explicitly
    
    invS = solve(S)
    ldS = log(det(S));  
    
    
    invSg = -invS%*%Sg%*%invS;
    ldSg = sum(diag(invS%*%Sg));
  }
  
  
  # update lambda
  if (is.na(eta.fixed)){
    eta = my$eta - alpha*Jg/p;
  } else {eta = my$eta}
  eta = min(eta,eta.max);
  eta = max(eta,eta.min);
  
  # return updated parameters
  
  my <- list(
    mu = mu, S = S, invS = invS, ldS = ldS, J = J,
    mug=mug, Sg=Sg, invSg=invSg, ldSg=ldSg, Jg=Jg,
    eta=eta
  )
  
  
  
  
  return(my)
}
