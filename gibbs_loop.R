
gibbs_loop <- function(n,ngibbs,b0,b1,signma2,plot=F){
  
  library(coda) 
  library(mvtnorm)
  
  beta <- matrix(c(b0,b1),2,1)		## put “true” regression parameters in a matrix
  
  x1 <- runif(n,0,20)
  x <- cbind(rep(1,n),x1)
  y <- matrix(rnorm(n,x%*%beta,sqrt(sigma2)),n,1)
  
  if(plot){plot(x1,y)
  abline(b0,b1,col=2,lwd=3)
  }
  
  #### specify priors
  bprior <- as.vector(c(0,0))
  vinvert <- solve(diag(1000,2))
  s1 <- 0.1
  s2 <- 0.1
  
  #### precompute frequently used quantities
  XX <- t(x) %*% x
  XY <- t(x) %*% y
  VbB <- vinvert %*% bprior
  
  #### storage for MCMC
  bgibbs <- matrix(0.0,nrow=ngibbs,ncol=2) 	## storage for beta
  sgibbs <- numeric(ngibbs)			## storage for sigma2
  
  #### initial conditions
  sg <- 50
  sinv <- 1/sg
  
  #### Gibbs Loop
  for(g in 1:ngibbs){
    
    ## sample regression parameters
    bigV    <- solve(sinv*XX + vinvert)  ## Covariance matrix
    littlev <- sinv*XY + VbB
    b = t(rmvnorm(1,bigV %*% littlev,bigV))   ## Vv is the mean vector
    
    ## sample variance
    u1 <- s1 + n/2
    u2 <- s2 + 0.5*crossprod(y-x%*%b)
    sinv <- rgamma(1,u1,u2)
    sg <- 1/sinv
    
    ## storage
    bgibbs[g,] <- b  ## store the current value of beta vector
    sgibbs[g]  <- sg  ## store the current value of the variance
    
    # if(g %%100 == 0) print(g)  ##show how many steps have been performed
  }
  
  out <- list(bgibbs=bgibbs,sgibbs=sgibbs,x1=x1,y=y)
  return(out)
}