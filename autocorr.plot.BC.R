autocorr.plot.BC <- function (x, lag.max, auto.layout = TRUE, ask, ...) 
{
  if (missing(ask)) {
    ask <- if (is.R()) {
      dev.interactive()
    }
    else {
      interactive()
    }
  }
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) 
    oldpar <- par(mfrow = set.mfrow(Nchains = nchain(x), 
                                    Nparms = nvar(x)))
  if (!is.mcmc.list(x)) 
    x <- mcmc.list(as.mcmc(x))
  for (i in 1:nchain(x)) {
    xacf <- if (missing(lag.max)) 
      acf(as.ts.mcmc(x[[i]]), plot = FALSE)
    else acf(as.ts.mcmc(x[[i]]), lag.max = lag.max, plot = FALSE)
    for (j in 1:nvar(x)) {
      plot(xacf$lag[, j, j], xacf$acf[, j, j], type = "h", 
           ylab = "Autocorrelation", xlab = "Lag", ylim = c(-1, 
                                                            1), ...)
      abline(h=.02,col="blue")
      title(paste(varnames(x)[j], ifelse(is.null(chanames(x)), 
                                         "", ":"), chanames(x)[i], sep = ""))
      if (i == 1 && j == 1) 
        oldpar <- c(oldpar, par(ask = ask))
    }
  }
  invisible(x)
}

