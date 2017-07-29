Stream.Gaussian <- function(
  x.new=x.new, my=list(), 
  eta.fixed=NA, alpha=0.001, eta.min=0.001, eta.max=0.4, efficient=FALSE, scaleit=FALSE
  ){

  ## x.new must be a vector
  ## params must be a list

  p = length(x.new); 

  ## INITIALISATIONS if t=0
  if (length(my)==0){
  	

		
    ## primary initialisations
    v.0 <- 10
    mu <- rep(0,p) 
    S <- diag(v.0,p,p)
    Pi <- S
    invS <- diag(1/v.0,p,p)
    ldS <- p*log(v.0) 
    J <- 0

    # gradients
    mug <- rep(0,p)
    Pig <- matrix(0,p,p)
    Sg <- matrix(0,p,p)
    invSg <- matrix(0,p,p) 
    ldSg <- 0 
    Jg <- 0 

    ## if adaptive forgetting, lambda=lambda.max, else lambda=lambda.fixed 
    if (is.na(eta.fixed))
    {eta<-eta.min} 
    else {eta<-eta.fixed}
    
    my <- list(
      mu = mu, S = S, Pi = Pi, invS = invS, ldS = ldS, J = J,
      mug=mug, Sg=Sg, Pig=Pig, invSg=invSg, ldSg=ldSg, Jg=Jg,
      eta=eta
    )

  }
 

  ## COMPUTE GRADIENT 
  J = my$ldS + t(x.new - my$mu)%*%my$invS%*%(x.new - my$mu)
  Jg = my$ldSg -2*t(x.new-my$mu)%*%(my$invS)%*%(my$mug) + t(x.new-my$mu)%*%my$invSg%*%(x.new-my$mu)
  if (scaleit) {Jg = (1/J)*Jg}

  ## PARAMETER UPDATES
  mu = my$mu + my$eta*(x.new - my$mu)
  mug = (1-my$eta)*my$mug +(x.new-my$mu);   
  
  Pi = my$Pi + my$eta*(x.new%*%t(x.new) - my$Pi) 
  Pig = (1-my$eta)*my$Pig +(x.new%*%t(x.new)-my$Pi);

  S = Pi - mu%*%t(mu)
  Sg = Pig - mug%*%t(mu) - mu%*%t(mug)

  if (efficient){
  	### NOT CHECKED YET !!!!!!
    # auxillary quantities gamma, h and their gradients
    g = (n-1)^2/n + t(xn-mu)*my$invS*(xn-mu)
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
    eta = my$eta - alpha*Jg;
  } else {eta = my$eta}
  eta = min(eta,eta.max);
  eta = max(eta,eta.min);

  # return updated parameters

  my <- list(
    mu = mu, S = S, Pi = Pi, invS = invS, ldS = ldS, J = J,
    mug=mug, Sg=Sg, Pig=Pig, invSg=invSg, ldSg=ldSg, Jg=Jg,
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

gen.univ <- function(true.params,n){
	
	
	true.params$rho.h = true.params$sd.h^2/(1-true.params$phi.h^2)
	true.params$rho.mu = true.params$sd.mu^2/(1-true.params$phi.mu^2)

	htrue <- rep(0,n)
	mutrue <- rep(0,n)
	htrue[1] <- rnorm(1,mean=0,sd=sqrt(true.params$rho.h))
	mutrue[1] <- rnorm(1,mean=0,sd=sqrt(true.params$rho.mu))
	for (t in 2:n){
		htrue[t] <- true.params$mu.h + true.params$phi.h*(htrue[t-1]-true.params$mu.h)+
		  rnorm(1,mean=0,sd=true.params$sd.h)
		mutrue[t] <- true.params$phi.mu*mutrue[t-1]+rnorm(1,mean=0,sd=true.params$sd.mu)
	}
	vartrue <- exp(htrue)*true.params$sd.v^2
	ys <- exp(htrue/2)*rnorm(n,mean=mutrue,sd=true.params$sd.v)
	return(list(ys=as.matrix(ys),mutrue=as.matrix(mutrue),vartrue=as.matrix(vartrue)))
}

gen.multiv <- function(true.params,n){

	p <- true.params$p
	
	df <- true.params$df*p
	Phi <- matrix(NA,p,p)
	for (i in 1:p){
		for (j in 1:p){
			if (abs(j-i)==1){
				Phi[i,j] <- -0.1
			}
			if (j==i){
				Phi[i,j] <- true.params$phi.h
			}
		}
	}
	
	mutrue <- matrix(NA,n,p)
	vartrue <- list()
	speed <- rep(NA,n)
	SNR <- rep(NA,n)
	ys <- matrix(NA,n,p)

	mutrue[1,] <- rnorm(p,mean=0,sd=1)
	vartrue[[1]] <- rwishart(df,(1/df)*diag(p))$W
	ys[1,] <- chol(vartrue[[1]])%*% rnorm(p,0,1)
	
	vartrue.vec <- matrix(NA,n,p^2)
	vartrue.vec[1,] <- t(as.vector(vartrue[[1]]))
	
	speed[1] <- p
	SNR[1] <- p
	Sigma.init <- (1/df)*diag(p)
	traces <- rep(0,n)
	for (t in 2:n){
	
#		vartrue[[t]] <- rwishart(df,(1/df)*diag(p))$W + rwishart(df,(1/df)*0.99*vartrue[[t-1]])$W
		if (true.params$phi.h==1){
			vartrue[[t]] <- vartrue[[t-1]]
		} else {
			vartrue[[t]] <- rwishart(df,(1/df)*Sigma.init+(1/df)*Phi*vartrue[[t-1]]*t(Phi))$W	
		}
		
		speed[t] <- sum(diag(solve(vartrue[[t]])%*%vartrue[[t-1]]))
		SNR[t] <- sum(diag(solve(vartrue[[t]])%*%Phi))
		vartrue.vec[t,] <- t(as.vector(vartrue[[t]]))
	
		traces[t-1] <- sum(diag(vartrue[[t]]%*%solve(vartrue[[t-1]])))
	
		mutrue[t,] <- true.params$phi.mu*mutrue[t-1,]+rnorm(p,mean=0,sd=true.params$sd.mu)
		ys[t,] <- chol(vartrue[[t]])%*% rnorm(p,0,1)
	}
	speed.scalar <- mean(log(speed))
	return(list(ys=as.matrix(ys),mutrue=as.matrix(mutrue),
	  vartrue=vartrue,vartrue.vec=vartrue.vec,speed=speed,SNR=SNR))
}

my.ellipse <- function(out.gen,nincr=30,m=4,m2=6,
  xlims=c(-50,50),ylims=c(-50,50)){
		library(ellipse)	
		M <- m*m2
		
		cols <- gray(rev((1:m2)/(m2+1)))
		ncount <- 1
		par(mfrow=c(sqrt(m),sqrt(m)))
		for (j in 1:m){
			plot(ellipse(out.gen$vartrue[[ncount]]),
			xlim=xlims,ylim=ylims,type='l',col=cols[1],
			main=paste('from n =',ncount,'to',ncount+m2*nincr))

			
			for (k in 1:m2){
				ncount <- ncount + nincr
				lines(ellipse(out.gen$vartrue[[ncount]]),
				xlim=xlims,ylim=ylims,col=cols[k+1])
					
			}
			
		}

}

my.multiv.plot <- function(out,ind.by=100,ind.j=c(1,2)){
	thelims <- c(
	  min(min(out$ys)),max(max(out$ys))
	)

	par(mfrow=c(3,3))
	for (i in 1:3){
		for (j in 1:3){
			k <- 3*(i-1)+j
			ind <- ((k-1)*ind.by+1):(k*ind.by)
			plot(out$ys[ind,ind.j[1]],out$ys[ind,ind.j[2]],
			  ylim=thelims,xlim=thelims,
			  xlab=paste('Variable ',ind.j[1],sep=''),
			  ylab=paste('Variable ',ind.j[2],sep=''),
			  main=paste('Time: ',min(ind),' to ',max(ind),sep='')
			)
			lines(ellipse(out$vartrue[[k*ind.by-ind.by/2]][ind.j,ind.j]))
		}
	}
	
	
}

run.experiment <- function(true.params,experiment.params,seed=1){

	set.seed(seed)

	# range of fixed learning rates to consider
	etarange = seq(experiment.params$etamin,
		experiment.params$etamax,
		experiment.params$etastep)
	m = length(etarange)

	alpha=experiment.params$alpha/p
	n=experiment.params$n
	N=experiment.params$N
	p=experiment.params$p

	# initialise (adaptive learning rates, MSE for adaptive, MSE for fixed)
	eta.adapt <- matrix(0,n,N)
	lds.adapt <- matrix(0,n,N)
	KL.adapt <- matrix(0,n,N)
	KL.fixed <- matrix(0,n,m)


	scount <- 1

	for (seed in 1:N){
	
		cat('Seed: ',seed,'\n')
	
		if (experiment.params$p==1){
			out.gen <- gen.univ(true.params,n=experiment.params$n+1)	
		} else {
			out.gen <- gen.multiv(true.params,n=experiment.params$n+1)	
		}
		
		ys <- out.gen$ys
		mutrue <- out.gen$mutrue
		vartrue <- out.gen$vartrue
		
		# RUN ADAPTIVE ALGORITHM
		params <- list()
		for (t in 1:n){
			params <- Stream.Gaussian(ys[t,],my=params,alpha=alpha)	
			eta.adapt[t,seed] <- params$eta
			lds.adapt[t,seed] <- params$ldS
			KL.adapt[t,seed] <- KL(mutrue[t,],params$mu,vartrue[[t]],params$S)
		}

		# OPTIONALLY RUN A NUMBER OF FIXED RATE ALGORITHMS
		if (seed <= experiment.params$nofixed){

			cat('Running fixed learning rates for seed ',seed,'...\n')				

			for (count in 1:m){			
				cat('eta = ',etarange[count],'\n')
				params <- list()
				for (t in 1:n){
					params <- Stream.Gaussian(ys[t,],my=params,eta.fixed=etarange[count])	
					KL.fixed[t,count] <- 
						(1/scount)*KL(mutrue[t,],params$mu,vartrue[[t]],params$S) 
						+ (1-1/scount)*KL.fixed[t,count]
				}
			}	
			scount <- scount + 1
		}
	}
  
  print(scount)
	out <- list(KL.fixed=KL.fixed,KL.adapt=KL.adapt,eta.adapt=eta.adapt,lds.adapt=lds.adapt)
	attr(out,'truth') <- out.gen
	return(out)
}


post.process <- function(true.params,experiment.params,out){
	
tolog <- experiment.params$tolog
if (tolog){

	KL.fixed <- log(out$KL.fixed)
	KL.adapt <- log(out$KL.adapt)

} else {

	KL.fixed <- (out$KL.fixed)
	KL.adapt <- (out$KL.adapt)

}

SNR <- round((true.params$sd.h/true.params$sd.v)^2,3)
speed <- round(true.params$sd.h^2,3)
if (experiment.params$p>1){
	speed <- round(mean(log(attr(out,'truth')$speed)),3)
	themain=paste(
	  'speed = ',speed,
	  ', alpha = ',alpha,sep=''
	)
} else {
	themain=paste(
	  'SNR = ',SNR,
	  ', speed = ',speed,
	  ', alpha = ',alpha,sep=''
	)
}
etarange = seq(experiment.params$etamin,experiment.params$etamax,experiment.params$etastep)
KL.fixed.all <- colMeans(KL.fixed[-c(1:n/5),])

eta.star <- etarange[which.min(KL.fixed.all)]
eta.best <- mean(mean(out$eta.adapt[-c(1:(n/5)),]))
KL.adapt.all <- rowMeans(KL.adapt)
KL.ensemble <- colMeans(KL.adapt[-c(1:(n/5)),])
KL.best <- mean(KL.ensemble)
KL.lb <- mean(KL.ensemble)-1.96*sqrt(var(KL.ensemble))
KL.ub <- mean(KL.ensemble)+1.96*sqrt(var(KL.ensemble))


thelims <- c(min(etarange),max(etarange))
themus <- rowMeans(out$eta.adapt)
stds <- sqrt(apply(out$eta.adapt,1,var))
theub <- themus + 1.96*stds
thelb <- themus - 1.96*stds

ylim <- c(
  min(min(KL.best,KL.lb),KL.fixed.all),
  max(max(KL.best,KL.ub),KL.fixed.all)
)

}
