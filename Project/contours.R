library(KernSmooth)

load("Project/mus.Rdata")

v1 = c(); v2 = c();
for(i in 1:5){
  v1 = c(v1, rep(i, 6-i))
  v2 = c(v2, (i+1):6)
}
  
par(mfcol = c(1,1))
k = length(out1.df[,1])

DEN.1 <- list()
for(N in 1:length(v1)){
  x <- out1.df[,v1[N]][1:k]
  y <- out1.df[,v2[N]][1:k]
  den <- bkde2D(cbind(x,y), bandwidth = c(bandwidth.nrd(x),bandwidth.nrd(y)))
  DEN.1 <- append(DEN.1, list(den))
}
DEN.2 <- list()
for(N in 1:length(v1)){
  x <- out2.df[,v1[N]][1:k]
  y <- out2.df[,v2[N]][1:k]
  den <- bkde2D(cbind(x,y), bandwidth = c(bandwidth.nrd(x),bandwidth.nrd(y)))
  DEN.2 <- append(DEN.2, list(den))
}
DEN.3 <- list()
for(N in 1:length(v1)){
  x <- out3.df[,v1[N]][1:k]
  y <- out3.df[,v2[N]][1:k]
  den <- bkde2D(cbind(x,y), bandwidth = c(bandwidth.nrd(x),bandwidth.nrd(y)))
  DEN.3 <- append(DEN.3, list(den))
}
DEN.4 <- list()
for(N in 1:length(v1)){
  x <- out4.df[,v1[N]][1:k]
  y <- out4.df[,v2[N]][1:k]
  den <- bkde2D(cbind(x,y), bandwidth = c(bandwidth.nrd(x),bandwidth.nrd(y)))
  DEN.4 <- append(DEN.4, list(den))
}

par(mfrow=c(2,2))
for(N in 1:15){
  z.1 <- DEN.1[[N]]$fhat; CI.1 <- quantile(z.1,c(.5)) 
  z.2 <- DEN.2[[N]]$fhat; CI.2 <- quantile(z.2,c(.5)) 
  z.3 <- DEN.3[[N]]$fhat; CI.3 <- quantile(z.3,c(.5)) 
  z.4 <- DEN.4[[N]]$fhat; CI.4 <- quantile(z.4,c(.5)) 

  xr <- range(c(DEN.1[[N]]$x1, DEN.3[[N]]$x1, DEN.2[[N]]$x1, DEN.4[[N]]$x1))
  yr <- range(c(DEN.1[[N]]$x2, DEN.3[[N]]$x2, DEN.2[[N]]$x2, DEN.4[[N]]$x2))
  
  pdfx.1 <- density(out1.df[,v1[N]])
  pdfx.2 <- density(out2.df[,v1[N]])
  pdfx.3 <- density(out3.df[,v1[N]])
  pdfx.4 <- density(out4.df[,v1[N]])
  pdfy.1 <- density(out1.df[,v2[N]])
  pdfy.2 <- density(out2.df[,v2[N]])
  pdfy.3 <- density(out3.df[,v2[N]])
  pdfy.4 <- density(out4.df[,v2[N]])
  
  contour(DEN.1[[N]]$x1,DEN.1[[N]]$x2,DEN.1[[N]]$fhat,levels=CI.1, labels=c(""), xlim=c(xr[1]-.1,xr[2]), ylim=c(yr[1]-.1,yr[2]), lty=2, lwd=2, col="red")
  contour(DEN.2[[N]]$x1,DEN.2[[N]]$x2,DEN.2[[N]]$fhat,levels=CI.3, labels=c(""), add=TRUE, lty=2, lwd=2, col="blue")
  contour(DEN.3[[N]]$x1,DEN.3[[N]]$x2,DEN.3[[N]]$fhat,levels=CI.2, labels=c(""), add=TRUE, lty=1, lwd=2, col="red")
  contour(DEN.4[[N]]$x1,DEN.4[[N]]$x2,DEN.4[[N]]$fhat,levels=CI.4, labels=c(""), add=TRUE, lty=1, lwd=2, col="blue")
  par(new=TRUE)
  plot(pdfx.1$x, pdfx.1$y, col="red", xlim=c(xr[1]-.1,xr[2]), ylim=c(0,15*max(pdfx.1$y)), type="l", lty=2,axes = FALSE, bty = "n", xlab = "", ylab = "")
  lines(pdfx.2$x, pdfx.2$y, col="blue", xlim=c(xr[1]-.1,xr[2]), ylim=c(0,10*max(pdfx.2$y)), lty=2)
  lines(pdfx.3$x, pdfx.3$y, col="red", xlim=c(xr[1]-.1,xr[2]), ylim=c(0,10*max(pdfx.3$y)))
  lines(pdfx.4$x, pdfx.4$y, col="blue", xlim=c(xr[1]-.1,xr[2]), ylim=c(0,10*max(pdfx.4$y)))
  par(new=TRUE)
  plot(pdfy.1$y, pdfy.1$x, col="red", xlim=c(0,xr[2]*300), ylim=c(yr[1]-.1,yr[2]), type="l", lty=2,axes = FALSE, bty = "n", xlab = "", ylab = "")
  lines(pdfy.2$y, pdfy.2$x, col="blue",xlim=c(0,xr[2]*300), ylim=c(yr[1]-.1,yr[2]), lty=2)
  lines(pdfy.3$y, pdfy.3$x, col="red", xlim=c(0,xr[2]*300), ylim=c(yr[1]-.1,yr[2]))
  lines(pdfy.4$y, pdfy.4$x, col="blue", xlim=c(0,xr[2]*300), ylim=c(yr[1]-.1,yr[2]))
  
  title(paste(colnames(out1.df[,v1])[N],"vs",colnames(out1.df[,v2])[N]))
}





