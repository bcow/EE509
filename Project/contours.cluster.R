#!/usr/bin/env Rscript

library(MASS)

load("/usr2/collab/ecowdery/EE509/Project/mus.Rdata")

v1 = c(); v2 = c();
for(i in 1:5){
  v1 = c(v1, rep(i, 6-i))
  v2 = c(v2, (i+1):6)
}

par(mfcol = c(1,1))
k = nrow(out1.df)

DEN.1 <- list()
for(N in 1:length(v1)){
  x <- out1.df[,v1[N]][1:k]
  y <- out1.df[,v2[N]][1:k]
  den <- kde2d(x, y, n = length(x))
  DEN.1 <- append(DEN.1, list(den))
}
DEN.2 <- list()
for(N in 1:length(v1)){
  x <- out2.df[,v1[N]][1:k]
  y <- out2.df[,v2[N]][1:k]
  den <- kde2d(x, y, n = length(x))
  DEN.2 <- append(DEN.2, list(den))
}
DEN.3 <- list()
for(N in 1:length(v1)){
  x <- out3.df[,v1[N]][1:k]
  y <- out3.df[,v2[N]][1:k]
  den <- kde2d(x, y, n = length(x))
  DEN.3 <- append(DEN.3, list(den))
}
DEN.4 <- list()
for(N in 1:length(v1)){
  x <- out4.df[,v1[N]][1:k]
  y <- out4.df[,v2[N]][1:k]
  den <- kde2d(x, y, n = length(x))
  DEN.4 <- append(DEN.4, list(den))
}

par(mfrow=c(1,1))
for(N in 1:15){
  z.1 <- DEN.1[[N]]$z; CI.1 <- quantile(z.1,c(.5))
  z.2 <- DEN.2[[N]]$z; CI.2 <- quantile(z.2,c(.5))
  z.3 <- DEN.3[[N]]$z; CI.3 <- quantile(z.3,c(.5))
  z.4 <- DEN.4[[N]]$z; CI.4 <- quantile(z.4,c(.5))
  
  xrange <- range(c(DEN.1[[N]]$x, DEN.3[[N]]$x, DEN.2[[N]]$x, DEN.4[[N]]$x))
  yrange <- range(c(DEN.1[[N]]$y, DEN.3[[N]]$y, DEN.2[[N]]$y, DEN.4[[N]]$y))
  
  contour(DEN.1[[N]],levels=CI.1, labels=c("97.5%"), xlim=xrange, ylim=yrange, lty=2, col="red")
  contour(DEN.3[[N]],levels=CI.3, labels=c("97.5%"), add=TRUE, lty=2, col="blue")
  contour(DEN.2[[N]],levels=CI.2, labels=c("97.5%"), add=TRUE, lty=1, col="red")
  contour(DEN.4[[N]],levels=CI.4, labels=c("97.5%"), add=TRUE, lty=1, col="blue")
}

save(DEN.1, DEN.2, DEN.3, DEN.4, file="/usr2/collab/ecowdery/EE509/Project/density.Rdata")