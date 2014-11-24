#########################
## Set up data
library(rjags)
glopnet <- read.csv("~/R Projects/EE509/Project/glopnet.csv")
glopnet[c("X")] <- NULL
glopnet[c("X.1")] <- NULL
glopnet[c("X.2")] <- NULL
glopnet[c("X.3")] <- NULL

gdata.na <- as.data.frame(cbind(
  glopnet$log.LL,
  glopnet$log.LMA,
  glopnet$log.Amass,
  glopnet$log.Nmass,
  glopnet$log.Pmass,
  glopnet$log.Rdmass
  ))
colnames(gdata.na)<-c("Log.LL","Log.LMA","Log.Amass","Log.Nmass","Log.Pmass","Log.Rmass")

gdata <- na.omit(gdata.na)

pairs(gdata.na)
pairs(gdata)
var(gdata.na,na.rm=TRUE)

############################
## MULTIVARIATE

MultModel = "
model{
  prec.Sigma~dwish(Vsig[,],n)
  Sigma[1:n,1:n] <- inverse(prec.Sigma[,])

  mu[1:n]~dmnorm(mu0[],Vmu)

  for(i in 1:N){
   for(j in 1:n){
     X[i,j]~dnorm(Y[i,j],10000000)
   }
  Y[i,1:n]~dmnorm(mu[],prec.Sigma[,])
  }
}"

# Without na's 
j.data <- gdata
N=dim(j.data)[1]; n=dim(j.data)[2]
data = list(Y=j.data, N=N, n=n, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
init = NULL
j.model   <- jags.model (file = textConnection(MultModel),data = data,inits = init,n.chains = 3)
update(j.model, n.iter=1000)
j.out   <- coda.samples (model = j.model,variable.names= c("mu"),n.iter = 10000)
out1 = j.out
out1.df <- as.data.frame(as.matrix(out1))

# With na's   ?????????????
# The problem is that a node in JAGS cannot be partially observed, whereas only the first element of X_mz_twin1[1,1:2] is observed.
# http://sourceforge.net/p/mcmc-jags/discussion/610037/thread/bad6e251/?limit=50
# Multivariate nodes cannot be partially observed, so if a node takes up two or more elements
#   of a node array, then the corresponding data values must be all present or all missing.

j.data <- gdata.na
N=dim(j.data)[1]; n=dim(j.data)[2]
data = list(X=j.data, N=N, n=n, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
init = NULL
j.model   <- jags.model (file = textConnection(MultModel), data = data, inits = init, n.chains = 3)
update(j.model, n.iter=1000)
j.out   <- coda.samples (model = j.model, variable.names= c("mu"), n.iter = 10000)
out2 = j.out
out2.df <- as.data.frame(as.matrix(out2))
 
#########
## SINGLE VARIATE

UnivModel = "
model{
  prec.sigma~dgamma(.001,.001)
  sigma <- 1/prec.sigma
  for(i in 1:n){mu[i]~dnorm(0,.001)}

  for(i in 1:N){
    for(j in 1:n){
      Y[i,j]~dnorm(mu[j],prec.sigma)
    }
  }
}"

# Without na's
j.data <- gdata
N=dim(j.data)[1]; n=dim(j.data)[2]
data = list(Y=j.data, N=N, n=n)
init = NULL
j.model   <- jags.model (file = textConnection(UnivModel),data = data,inits = init, n.chains = 3)
update(j.model, n.iter=1000)
j.out   <- coda.samples (model = j.model,variable.names= c("mu"), n.iter = 10000)
out3 = j.out
out3.df <- as.data.frame(as.matrix(out3))

# With na's
j.data <- gdata.na
N=dim(j.data)[1]; n=dim(j.data)[2]
data = list(Y=j.data, N=N, n=n)
init = NULL
j.model   <- jags.model (file = textConnection(UnivModel),data = data,inits = init, n.chains = 3)
update(j.model, n.iter=1000)
j.out   <- coda.samples (model = j.model,variable.names= c("mu"), n.iter = 10000)
out4 = j.out
out4.df <- as.data.frame(as.matrix(out4))

############################
# Comparison
rbind(
  colMeans(j.data),
  c(mean(out1.df$"mu[1]"),mean(out1.df$"mu[2]"),mean(out1.df$"mu[3]"),mean(out1.df$"mu[4]"),mean(out1.df$"mu[5]"),mean(out1.df$"mu[6]")),
  c(mean(out2.df$"mu[1]"),mean(out2.df$"mu[2]"),mean(out2.df$"mu[3]"),mean(out2.df$"mu[4]"),mean(out2.df$"mu[5]"),mean(out2.df$"mu[6]")),
  c(mean(out3.df$"mu[1]"),mean(out3.df$"mu[2]"),mean(out3.df$"mu[3]"),mean(out3.df$"mu[4]"),mean(out3.df$"mu[5]"),mean(out3.df$"mu[6]")),
  c(mean(out4.df$"mu[1]"),mean(out4.df$"mu[2]"),mean(out4.df$"mu[3]"),mean(out4.df$"mu[4]"),mean(out4.df$"mu[5]"),mean(out4.df$"mu[6]"))
)

par(mfrow=c(1,1))
par(mfrow = c(3,2))



ranges <-  apply(cbind(out1.df$"mu[1]",out2.df$"mu[1]",out3.df$"mu[1]",out4.df$"mu[1]"), 2,
                 function(x) { dens <- density(x); c(range(dens$x), range(dens$y)) })
plot(density(out1.df$"mu[1]"), col="red", 
     xlim = range(ranges[1:2, ]), ylim = range(ranges[3:4, ]),
     main=paste(names(gdata)[1]))
lines(density(out2.df$"mu[1]"), col="blue")
lines(density(out3.df$"mu[1]"), col="red", lty=2)
lines(density(out4.df$"mu[1]"), col="blue", lty=2)

ranges <-  apply(cbind(out1.df$"mu[2]",out2.df$"mu[2]]",out3.df$"mu[2]",out4.df$"mu[2]"), 2,
                 function(x) { dens <- density(x); c(range(dens$x), range(dens$y)) })
plot(density(out1.df$"mu[2]"), col="red", 
     xlim = range(ranges[1:2, ]), ylim = range(ranges[3:4, ]),
     main=paste(names(gdata)[2]))
lines(density(out2.df$"mu[2]"), col="blue")
lines(density(out3.df$"mu[2]"), col="red", lty=2)
lines(density(out4.df$"mu[2]"), col="blue", lty=2)

ranges <-  apply(cbind(out1.df$"mu[3]",out2.df$"mu[3]",out3.df$"mu[3]",out4.df$"mu[3]"), 2,
                 function(x) { dens <- density(x); c(range(dens$x), range(dens$y)) })
plot(density(out1.df$"mu[3]"), col="red", 
     xlim = range(ranges[1:2, ]), ylim = range(ranges[3:4, ]),
     main=paste(names(gdata)[3]))
lines(density(out2.df$"mu[3]"), col="blue")
lines(density(out3.df$"mu[3]"), col="red", lty=2)
lines(density(out4.df$"mu[3]"), col="blue", lty=2)

ranges <-  apply(cbind(out1.df$"mu[4]",out2.df$"mu[4]",out3.df$"mu[4]",out4.df$"mu[4]"), 2,
                 function(x) { dens <- density(x); c(range(dens$x), range(dens$y)) })
plot(density(out1.df$"mu[4]"), col="red", 
     xlim = range(ranges[1:2, ]), ylim = range(ranges[3:4, ]),
     main=paste(names(gdata)[4]))
lines(density(out2.df$"mu[4]"), col="blue")
lines(density(out3.df$"mu[4]"), col="red", lty=2)
lines(density(out4.df$"mu[4]"), col="blue", lty=2)

ranges <-  apply(cbind(out1.df$"mu[5]",out2.df$"mu[5]",out3.df$"mu[5]",out4.df$"mu[5]"), 2,
                 function(x) { dens <- density(x); c(range(dens$x), range(dens$y)) })
plot(density(out1.df$"mu[5]"), col="red", 
     xlim = range(ranges[1:2, ]), ylim = range(ranges[3:4, ]),
     main=paste(names(gdata)[5]))
lines(density(out2.df$"mu[5]"), col="blue")
lines(density(out3.df$"mu[5]"), col="red", lty=2)
lines(density(out4.df$"mu[5]"), col="blue", lty=2)

ranges <-  apply(cbind(out1.df$"mu[6]",out2.df$"mu[6]",out3.df$"mu[6]",out4.df$"mu[6]"), 2,
                 function(x) { dens <- density(x); c(range(dens$x), range(dens$y)) })
plot(density(out1.df$"mu[6]"), col="red", 
     xlim = range(ranges[1:2, ]), ylim = range(ranges[3:4, ]),
     main=paste(names(gdata)[6]))
lines(density(out2.df$"mu[6]"), col="blue")
lines(density(out3.df$"mu[6]"), col="red", lty=2)
lines(density(out4.df$"mu[6]"), col="blue", lty=2)

#####################
## Comparison Tables

# 

















# 
# 
# #############
# ## With missing values
# 
# j.data <- gdata.na[1:50,1:5]
# N=dim(j.data)[1]; n=dim(j.data)[2]
# names <- names(j.data)
# 
# na.idx <- list() 
# for(name in names){
#   na.idx <- append(na.idx, list(which(is.na(j.data[name]))))
# }
# names(na.idx) <- paste0("NA",1:n)
# #
# 
# data = append(list(Y=j.data, N=N, n=n, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n)), na.idx)
# init = NULL
# 
# 
# 
# # for(i in 1:N){
# #   for(j in 1:n){
# #     Y[i,j]~dnorm(0,.001)
# #   }
# # }
# 
# 
# # for(i in 1:length(NA1)){Y[NA1[i],1]~dnorm(0,.001)}
# # for(i in 1:length(NA2)){Y[NA2[i],2]~dnorm(0,.001)}
# # for(i in 1:length(NA3)){Y[NA3[i],3]~dnorm(0,.001)}
# # for(i in 1:length(NA4)){Y[NA4[i],4]~dnorm(0,.001)}
# # for(i in 1:length(NA5)){Y[NA5[i],5]~dnorm(0,.001)}
# # for(i in 1:length(NA6)){Y[NA6[i],6]~dnorm(0,.001)}
# 
# 
# 
# 
# 
# FitNorm = "
# model{
# 
# ym~dnorm(0,.001)
# tau~dgamma(.001,.001)
# 
# for(i in 1:length(NA1)){Y[NA1[i],1]~dnorm(0,.001)}
# #for(i in 1:length(NA2)){Y[NA2[i],2]~dnorm(0,.001)}
# for(i in 1:length(NA3)){Y[NA3[i],3]~dnorm(0,.001)}
# for(i in 1:length(NA4)){Y[NA4[i],4]~dnorm(0,.001)}
# for(i in 1:length(NA5)){Y[NA5[i],5]~dnorm(0,.001)}
# 
#   prec.Sigma~dwish(Vsig[,],n)
#   Sigma[1:n,1:n] <- inverse(prec.Sigma[,])
#   
#   mu[1:n]~dmnorm(mu0[],Vmu)
#   
#   for(i in 1:N){
#     Y[i,1:n]~dmnorm(mu[],prec.Sigma[,])
#   }
# }"
# 
# j.model   <- jags.model (file = textConnection(FitNorm),
#                          data = data,
#                          inits = init,
#                          n.chains = 3)
# update(j.model, n.iter=1000)
# j.out   <- coda.samples (model = j.model,
#                          variable.names= c("mu", "Sigma"),
#                          n.iter = 10000)
# 
# out1 = j.out
# summary(out1)
# plot(out1)
# autocorr.plot(out1)
# gelman.plot(out1) 
# 
# cumuplot(out1)
# 
# out1.df <- as.data.frame(as.matrix(out1))
# m <- matrix(out1.df,nrow = 72)
# m.array<-array(out1.df,dim=c(72,1,30000))
# 
# rbind(
#   colMeans(j.data),
#   c(mean(out1.df$"mu[1]"),mean(out1.df$"mu[2]"),mean(out1.df$"mu[3]"),mean(out1.df$"mu[4]"),mean(out1.df$"mu[5]"),mean(out1.df$"mu[6]"))
# )
# 
# #############
# ## With only one missing value
# 
# j.data <- gdata.na
# N=dim(j.data)[1]; n=dim(j.data)[2]
# data = list(X=j.data, N=N, n=n, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
# init = NULL
# 
# FitNorm = "
# model{
# prec.Sigma~dwish(Vsig[,],n)
# Sigma[1:n,1:n] <- inverse(prec.Sigma[,])
# 
# mu[1:n]~dmnorm(mu0[],Vmu)
# 
# for(i in 1:N){
#   for(j in 1:n){
#     X[i,j]~dnorm(Y[i,j],10000000)
#   }
# Y[i,1:n]~dmnorm(mu[],prec.Sigma[,])
# }
# }"
# 
# j.model   <- jags.model (file = textConnection(FitNorm), data = data, inits = init, n.chains = 3)
# update(j.model, n.iter=1000)
# j.out   <- coda.samples (model = j.model, variable.names= c("mu", "Sigma"), n.iter = 10000)
# 
# out4 = j.out
# # summary(out1)
# plot(out4)
# # autocorr.plot(out1)
# # gelman.plot(out1) 
# # cumuplot(out1)
# 
# mu[1]       1.235368 0.030257 1.747e-04      1.736e-04
# mu[2]       2.203356 0.028255 1.631e-04      1.627e-04
# mu[3]       1.772627 0.029422 1.699e-04      1.680e-04
# mu[4]       0.138907 0.026078 1.506e-04      1.491e-04
# mu[5]      -1.257074 0.034276 1.979e-04      1.979e-04
# mu[6]       0.915528 0.033002 1.905e-04      1.905e-04
# 
# mu[1]       1.235044 0.030336 1.751e-04      1.751e-04
# mu[2]       2.203494 0.028297 1.634e-04      1.620e-04
# mu[3]       1.772793 0.029473 1.702e-04      1.693e-04
# mu[4]       0.138782 0.025916 1.496e-04      1.479e-04
# mu[5]      -1.257108 0.034360 1.984e-04      1.982e-04
# mu[6]       0.915613 0.033049 1.908e-04      1.891e-04
# 
# out1.df <- as.data.frame(as.matrix(out1))
# rbind(
#   colMeans(j.data),
#   c(mean(out1.df$"mu[1]"),mean(out1.df$"mu[2]"),mean(out1.df$"mu[3]"),mean(out1.df$"mu[4]"),mean(out1.df$"mu[5]"),mean(out1.df$"mu[6]"))
# )

# 
# ############
# ## Regression
# 
# library(lmodel2)
# 
# fit <- lmodel2(gdata$Log.Nmas s~ gdata$Log.Pmass)
# plot(fit)
# 
# var(gdata$Log.Pmass, gdata$Log.Nmass, na.rm=TRUE)
# 
# 
# FitNorm = "
# model{
# for(i in 1:3){b[i]~dnorm(0,.001)}
# 
# prec.sigma2~dgamma(.001,.001)
# sigma2 <- 1/prec.sigma2
# 
# prec.Sigma~dwish(V[,],3)
# Sigma[1:3,1:3] <- inverse(prec.Sigma[,])
# 
# for(i in 1:n){
# beta[i,1:3]~dmnorm(b[],prec.Sigma[,])
# mu[i] <- beta[i,1] + beta[i,2]*N[i] + beta[i,3]*P[i]
# LMA[i]~dnorm(mu[i],prec.sigma2)
# }
# }"
# 
# V <- diag(3)
# data = list(LMA=gdata$Log.LMA, N=gdata$Log.Nmass, P=gdata$Log.Pmass, n=length(gdata$Log.LMA), V=V)
# init = NULL
# 
# # init.cond <- list()
# # init.cond1[[1]] = list(b=c(1,1,1), prec.sigma2 = 1, prec.Sigma = structure(.Data=c(1,0,0,0,1,0,0,0,1),.Dim=c(3,3)))
# # init.cond1[[2]] = 
# # init.cond1[[3]] = 
# # init = init.cond1
# 
# 
# 
# j.model   <- jags.model (file = textConnection(FitNorm),
#                          data = data,
#                          inits = init,
#                          n.chains = 3)
# update(j.model, n.iter=1000)
# j.out   <- coda.samples (model = j.model,
#                          variable.names = c("sigma2","b","Sigma"),
#                          n.iter = 200000,
#                          thin = 40)
# 
# out1 = j.out
# summary(out1)
# plot(out1)
# autocorr.plot(out1)