# Looking at Sigma

library(MASS)

par(mfrow = c(6,6))
Sigma <- as.data.frame(matrix(NA,6,6))
for(i in 1:6){
  for(j in 1:6){
    Sigma[i,j] <- mean(out4.df[,i+(j-1)*6])
    plot(density(abs(out4.df[,i+(j-1)*6])),main=colnames(out4.df)[i+(j-1)*6], xlim=c(abs(var(gdata)[i,j])-.01, max(abs(out4.df[,i+(j-1)*6]))))
    abline(v=abs(Sigma[i,j]),col=3, lwd=2)
    abline(v=abs(var(gdata)[i,j]), lty=2, col=2, lwd=2)
  }
}
colnames(Sigma)<- rownames(Sigma) <- c("Log.LL","Log.LMA","Log.Amass","Log.Nmass","Log.Pmass","Log.Rmass")

print(Sigma)
print(var(gdata))
print(abs(Sigma) >= abs(var(gdata)))


pairs(out4.df[,37:42], panel=function(x,y){
  points(x,y)
  fit <- lm(y~x)
  p <- pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],summary(fit)$fstatistic[3], lower.tail = F)
  if(p < .01){abline(fit, col='red',lwd=2)}
  legend("top", legend=sprintf("R2 = %.2f",summary(fit)$r.squared), text.col="blue")
})

pairs(gdata, panel=function(x,y){
#  points(x,y)
#   fit <- lm(y~x)
#   p <- pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],summary(fit)$fstatistic[3], lower.tail = F)
#   if(p < .01){abline(fit, col='red',lwd=2)}
#   par(new=TRUE)
  den <- (kde2d(x, y, n = length(x)))
  z <- den$z
  CI <- quantile(z,c(.5,.95)) 
  contour(den, col = "red", levels=CI ,add = TRUE) 
  #legend("top", legend=sprintf("R2 = %.2f",summary(fit)$r.squared), text.col="blue")
})

pairs(out4.df[,37:42], panel=function(x,y){
  # points(x,y)
  #   fit <- lm(y~x)
  #   p <- pf(summary(fit)$fstatistic[1],summary(fit)$fstatistic[2],summary(fit)$fstatistic[3], lower.tail = F)
  #   if(p < .01){abline(fit, col='red',lwd=2)}
  #   par(new=TRUE)
  den <- (kde2d(x, y, n = length(x)))
  z <- den$z
  CI <- quantile(z,c(.5,.95)) 
  contour(den, col = "red", levels=CI ,add = TRUE) 
  #legend("top", legend=sprintf("R2 = %.2f",summary(fit)$r.squared), text.col="blue")
})
print((cov2cor(as.matrix(Sigma)))^2)
print((cor(gdata))^2)
print(abs((cov2cor(as.matrix(Sigma)))^2) >= abs((cor(gdata))^2))
