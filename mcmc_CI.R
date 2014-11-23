mcmc_CI <- function(ngibbs,bgibbs,sgibbs,beg,thin,x1,y,b0,b1,plot=F){
  
  xpred <- 0:20      			## sequence of x values we're going to
  npred <- length(xpred)				##      make predictions for
  ypred <- matrix(0.0,nrow=ngibbs,ncol=npred)	## storage for predictive interval
  ycred <- matrix(0.0,nrow=ngibbs,ncol=npred)	## storage for credible interval
  
  for(g in seq(from=beg,to=ngibbs,by=thin)){
    Ey <- bgibbs[g,1] + bgibbs[g,2] * xpred
    ycred[g,] <- Ey
    ypred[g,] <- rnorm(npred,Ey,sqrt(sgibbs[g]))
  }
  
  ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))  ## credible interval and median
  pi <- apply(ypred,2,quantile,c(0.025,0.975))  	## prediction interval
  
  if(plot){
    plot(x1,y,cex=0.5,xlim=c(0,20),ylim=c(0,50))
    lines(xpred,ci[1,],col=3,lty=2)	## lower CI
    lines(xpred,ci[2,],col=3,lwd=2)	## median
    lines(xpred,ci[3,],col=3,lty=2)	## upper CI
    lines(xpred,pi[1,],col=4,lty=2)	## lower PI
    lines(xpred,pi[2,],col=4,lty=2)	## upper PI
    abline(b0,b1)				## true model
  }
}

