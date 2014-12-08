library("rjags")
library(coda)
library(mvtnorm)
library(nnet)

setwd("/Users/elizabethcowdery/R Projects/EE509/Project")
glopnet <- read.csv("glopnet.csv")
glopnet[c("X")] <- NULL
glopnet[c("X.1")] <- NULL
glopnet[c("X.2")] <- NULL
glopnet[c("X.3")] <- NULL

DSI <- as.numeric(glopnet$Dataset)

gdata.na <- as.data.frame(cbind(
  glopnet$log.LL,
  glopnet$log.LMA,
  glopnet$log.Amass,
  glopnet$log.Nmass,
  glopnet$log.Pmass,
  glopnet$log.Rdmass 
  #   glopnet$Species, 
  #   glopnet$BIOME,
  #   glopnet$Dataset
))
colnames(gdata.na)<-c("Log.LL","Log.LMA","Log.Amass","Log.Nmass","Log.Pmass","Log.Rmass")
gdata <- na.omit(gdata.na)

MultModel = "
model{
  prec.Sigma~dwish(Vsig[,],n)
  Sigma[1:n,1:n] <- inverse(prec.Sigma[,])
  
  mu[1:n]~dmnorm(mu0[],Vmu)
  
  for(i in 1:N){
    Y[i,1:n]~dmnorm(mu[],prec.Sigma[,])
    for(j in 1:n){
      X[i,j]~dnorm(Y[i,j],10000000)
    }
  }
}"

# Without na's 
j.data <- gdata
N=dim(j.data)[1]; n=dim(j.data)[2]
data = list(Y=j.data, N=N, n=n, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
init = NULL
j.model   <- jags.model (file = textConnection(MultModel),data = data,inits = init,n.chains = 3)
# update(j.model, n.iter=1000)
# j.out   <- coda.samples (model = j.model,variable.names= c("mu"),n.iter = 10000)
# out2 = j.out


# With na's 

j.data <- gdata.na
N=dim(j.data)[1]; n=dim(j.data)[2]
data = list(X=j.data, N=N, n=n, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
init = NULL
j.model   <- jags.model (file = textConnection(MultModel), data = data, inits = init, n.chains = 3)
# update(j.model, n.iter=1000)
# j.out   <- coda.samples (model = j.model, variable.names= c("mu"), n.iter = 10000)
# out4 = j.out


RandMultModel = "
model{
  prec.Sigma~dwish(Vsig[,],n)
  Sigma[1:n,1:n] <- inverse(prec.Sigma[,])

  for(i in 1:nds){alpha.ds[i]~dnorm(0,tau.d)}
  tau.d~dgamma(.001,.001)

  for(i in 1:N){
    for(j in 1:n){
      a[i,j] <- alpha.ds[DSI[i]]
    }
  }

  beta[1:n]~dmnorm(mu0[],Vmu)

  for(i in 1:N){mu[i,1:n] <- beta[1:n]+a[i,1:n]}


  for(i in 1:N){
    Y[i,1:n]~dmnorm(mu[i,1:n],prec.Sigma[,])
    for(j in 1:n){
      X[i,j]~dnorm(Y[i,j],10000000)
    }
  }
}"

# Without na's 
j.data <- gdata.na
N=dim(j.data)[1]; n=dim(j.data)[2]; nds = length(unique(DSI))
data = list(Y=j.data, N=N, n=n, nds=nds, DSI=DSI, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
init = NULL
j.model   <- jags.model (file = textConnection(RandMultModel),data = data,inits = init,n.chains = 3)
# update(j.model, n.iter=1000)
# j.out   <- coda.samples (model = j.model,variable.names= c("mu"),n.iter = 10000)
# out = j.out


# With na's 

j.data <- gdata.na
DSI <- class.ind(j.data$Dataset)
j.data$Dataset <- NULL
N=dim(j.data)[1]; n=dim(j.data)[2]; nds = dim(DSI)[2]
data = list(Y=j.data, N=N, n=n, nds=nds, DSI=DSI, Vsig = diag(n), mu0 = rep(0,n), Vmu = diag(.001,n))
init = NULL
j.model   <- jags.model (file = textConnection(RandMultModel), data = data, inits = init, n.chains = 3)

update(j.model, n.iter=1000)
j.out   <- coda.samples (model = j.model, variable.names= c("mu"), n.iter = 10000)
out4 = j.out

out2.df <- as.data.frame(as.matrix(out2))
out4.df <- as.data.frame(as.matrix(out4))

colnames(out2.df)<-colnames(out4.df)<-colnames(gdata)











