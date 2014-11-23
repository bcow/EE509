mcmc_diag_MM <- function(out,beg,thin,plot=F){
  
  ## convert to MCMC object
  mcmc <- mcmc(cbind(
    out$bgibbs[seq(from=beg,to=ngibbs,by=thin),],
    out$sgibbs[seq(from=beg,to=ngibbs,by=thin)],
    out$tgibbs[seq(from=beg,to=ngibbs,by=thin)]
    )) 
  if(plot){
    print("var 1 = b0, var2 = b1, var3 = variance, v4 = theta")
    plot(mcmc) ## mcmc history and density plot
    
    autocorr.plot(mcmc)  	## autocorrelation
    cumuplot(mcmc)		## quantile plot
    
    print("var 1 = b0, var2 = b1, var3 = variance, v4 = theta")
    print(1-rejectionRate(mcmc))	## acceptance rate
    print(summary(mcmc))		## summary table
    
    par(mfrow = c(1,1))
    plot(out$bgibbs[,1],out$bgibbs[,2],xlab="b0",ylab="b1",main="b0 vs b1")	## pairs plot to evaluate parameter correlation
  }
  return(summary(mcmc))
}