
gibbs_loop_MM <- function(L,ngibbs,b0,b1,signma2,plot=F){
  
  library(coda) 
  library(mvtnorm)
  
  beta <- matrix(c(b0,b1),2,1)  	## put “true” regression parameters in a matrix
  
  #### specify priors
  bprior <- as.vector(c(0,0))
  vinvert <- solve(diag(1000,2))
  s1 <- 0.1
  s2 <- 0.1
  
  ### beta prior for theta
  a1 = 1.91
  a2 = 10.17
  
  Vb <- vinvert %*% bprior
  
  #### storage for MCMC
  bgibbs <- matrix(0.0,nrow=ngibbs,ncol=2) 	## storage for beta
  sgibbs <- numeric(ngibbs)			## storage for sigma2
  tgibbs <- rep(theta,ngibbs)   ## storage for theta
  

  #### initial conditions
  sg <- 50
  sinv = 1/sg
  n <- length(L)
  z <- L/(L+theta)
  X <- cbind(rep(1,n),z)

  Y <- matrix(rnorm(n,X%*%beta,sqrt(sigma2)),n,1)
  
  if(plot){
    plot(L,Y)
    lines(l,Mic.Men(beta,l,theta),col="red")
    points(L, grow, col="blue")
    abline(b0,b1,col=2,lwd=3)
  }
  
  
  #### Gibbs Loop
  for(g in 1:ngibbs){
    
    ## sample regression parameters
    bigV    <- solve(sinv*crossprod(X) + vinvert)
    littlev <- sinv*crossprod(X,grow) + Vb
    b <- t(rmvnorm(1,bigV %*% littlev,bigV))
    
    ## sample variance
    u1 <- s1 + n/2
    u2 <- s2 + 0.5*crossprod(grow-X%*%b)
    sinv <- rgamma(1,u1,u2)
    sg <- 1/sinv
    
    ##theta
    tnew <- rtnorm(1,theta,JumpSD)        ##propose new theta
    znew <- L/(L+tnew)                    ## calculate new z
    Xnew <- cbind(rep(1,n),znew)              ## calculate new X
    anum <- dmvnorm(grow,Xnew%*%b,diag(sg,n),log=TRUE) +  ##likelihood
      dbeta(tnew,a1,a2,log=TRUE)          ##prior
    jnum <- dtnorm(tnew,theta,JumpSD)             ##jump
    aden <- dmvnorm(grow,X%*%b,diag(sg,n),log=TRUE) + ##likelihood
      dbeta(theta,a1,a2,log=TRUE)           ##prior
    jden <- dtnorm(theta,tnew,JumpSD)             ##jump
    a <- exp((anum-jnum)-(aden-jden))         ## acceptance criteria
    if(a > runif(1)){                 ## accept with probability a
      theta <- tnew                       ## update theta if step accepted
      X <- Xnew                       ## update X if step accepted
    }
    
    ## storage
    bgibbs[g,] <- b  ## store the current value of beta vector
    sgibbs[g]  <- sg  ## store the current value of the variance
    tgibbs[g] <- tnew  ## store the current value of  theta
    
    # if(g %%100 == 0) print(g)  ##show how many steps have been performed
  }
  
  out <- list(bgibbs=bgibbs,sgibbs=sgibbs, tgibbs=tgibbs,X=X,Y=Y)
  return(out)
}