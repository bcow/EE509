

# Fit the overall mean and standard deviation, reporting summary statistics for both

# FitNorm = "
# model {
#   mu ~ dnorm(0,0.001)
#   prec ~ dgamma(0.1,0.1)
#     for(t in 1:nt){
#       for(i in 1:nrep){
#       x[t,i] ~ dnorm(mu,prec)
#     }
#   }
# }"
# 
# data = list(x=mos,nt=10,nrep=5)
# init = NULL
# 
# j.model   <- jags.model (file = textConnection(FitNorm),
#                          data = data,
#                          inits = init,
#                          n.chains = 3)
# b1   <- coda.samples (model = j.model,
#                       variable.names = c("mu","prec","x"),
#                       n.iter = 2000)
# bmcmc = b1
#  
# autocorr.plot(bmcmc)
# cumuplot(bmcmc)
# gelman.plot(bmcmc) 
# plot(bmcmc) 
# summary(bmcmc) 
# 
# mu = as.data.frame(as.matrix(bmcmc))$mu
# quantile(mu,c(0.025,0.5,0.975))


###################

library(rjags)

mosq <- read.csv(file = "data//Mosquito.csv")

ymax <- max(mosq$rep1,mosq$rep2,mosq$rep3,mosq$rep4,mosq$rep5)
ymin <- min(mosq$rep1,mosq$rep2,mosq$rep3,mosq$rep4,mosq$rep5)

# Plot mosquito abundance as a function of time in a way that distinguishes the reps 
plot(mosq$time,mosq$rep1, type="l", col=1, ylim=c(ymin, ymax))
lines(mosq$time,mosq$rep2, col=2)
lines(mosq$time,mosq$rep3, col=3)
lines(mosq$time,mosq$rep4, col=4)
lines(mosq$time,mosq$rep5, col=5)


mos <- cbind(mosq$rep1,mosq$rep2,mosq$rep3,mosq$rep4,mosq$rep5)
  
  
data = list(x=mos,nt=10,nrep=5)
init.cond1 <- list()
init.cond1[[1]] = list(mu=7)
init.cond1[[2]] = list(mu=7.5)
init.cond1[[3]] = list(mu=8)
init = init.cond1

FitNorm = "
model{
  mu ~ dnorm(0,0.001)       ## prior mean
  sigma ~ dgamma(0.001,0.001)   ## prior residual precision
  tau.t ~ dgamma(0.001,0.001)   ## prior year-effect precision
  
  for(t in 1:nt){           ## loop over years
    alpha.t[t] ~ dnorm(0,tau.t)     ## random year effect
    Ex[t] <- mu + alpha.t[t]        ## process model (does not vary with rep i)
    for(i in 1:nrep){           ## loop over reps
        x[t,i] ~ dnorm(Ex[t],sigma) ## data model
    }
    px[t] ~ dnorm(Ex[t],sigma)  ## predictive interval
  }
}
"

j.model   <- jags.model (file = textConnection(FitNorm),
                         data = data,
                         inits = init,
                         n.chains = 3)
b1   <- coda.samples (model = j.model,
                      variable.names = c("Ex","px","mu","sigma","tau.t","alpha.t"),
                      n.iter = 130000,
                      burnin=2000,
                      thin=25)
bmcmc1 = b1

autocorr.plot(bmcmc1)
cumuplot(bmcmc1)
gelman.plot(bmcmc1) 
plot(bmcmc1)  
summary(bmcmc1) 

bmcmc1.df <- as.data.frame(as.matrix(bmcmc1))

quants <- rbind(
quantile(bmcmc1.df$'Ex[1]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[2]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[3]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[4]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[5]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[6]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[7]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[8]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[9]',c(0.025,0.5,0.975)),
quantile(bmcmc1.df$'Ex[10]',c(0.025,0.5,0.975)))

preds <- rbind(
  quantile(bmcmc1.df$'px[1]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[2]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[3]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[4]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[5]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[6]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[7]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[8]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[9]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[10]',c(0.025,0.5,0.975)))


ymin=min(preds[,1]);ymax=max(preds[,3])
plot(mosq$time,mosq$rep1, type="l", col=8, ylim=c(ymin, ymax))
lines(mosq$time,mosq$rep2, col=8)
lines(mosq$time,mosq$rep3, col=8)
lines(mosq$time,mosq$rep4, col=8)
lines(mosq$time,mosq$rep5, col=8)
lines(mosq$time,quants[,1],lty=2, col=2, lwd=2)
lines(mosq$time,quants[,2],col=2, lwd=2)
lines(mosq$time,quants[,3],lty=2,col=2, lwd=2)
lines(mosq$time,preds[,1],lty=2, col=3, lwd=2)
lines(mosq$time,preds[,3],lty=2, col=3, lwd=2)

######################

pre.means <- c(
  mean(bmcmc1.df$'Ex[1]'),
  mean(bmcmc1.df$'Ex[2]'),
  mean(bmcmc1.df$'Ex[3]'),
  mean(bmcmc1.df$'Ex[4]'),
  mean(bmcmc1.df$'Ex[5]'),
  mean(bmcmc1.df$'Ex[6]'),
  mean(bmcmc1.df$'Ex[7]'),
  mean(bmcmc1.df$'Ex[8]'),
  mean(bmcmc1.df$'Ex[9]'),
  mean(bmcmc1.df$'Ex[10]'))

obs.means <- rowMeans(mos)

mm <- mean(mos)

ss.tot <- sum((obs.means-mm)^2)
ss.reg <- sum((pre.means-mm)^2)
ss.res <- sum((obs.means-pre.means)^2)

rsq <- ss.reg/ss.tot


########

alpha.t <- rowMeans(rbind(
  bmcmc1.df$'alpha.t[1]',
  bmcmc1.df$'alpha.t[2]',
  bmcmc1.df$'alpha.t[3]',
  bmcmc1.df$'alpha.t[4]',
  bmcmc1.df$'alpha.t[5]',
  bmcmc1.df$'alpha.t[6]',
  bmcmc1.df$'alpha.t[7]',
  bmcmc1.df$'alpha.t[8]',
  bmcmc1.df$'alpha.t[9]',
  bmcmc1.df$'alpha.t[10]'
             ))

plot(rowMeans(alpha.t))
met <- read.csv("data//met.csv")

time <- 1995:2004
idx <- which(met$year == time)
par(mfrow=c(4,1))
plot(time,alpha.t, type="o",lwd=2); par(new=TRUE)
plot(time, met$precip[idx], type="o",col=2)
plot(time,alpha.t, type="o",lwd=2); par(new=TRUE)
plot(time, met$MAT[idx], type="o",col=3)
plot(time,alpha.t, type="o",lwd=2); par(new=TRUE)
plot(time, met$RH[idx], type="o",col=4)


data = list(x=mos,precip=met$precip,nt=10,nrep=5)
init.cond1 <- list()
init.cond1[[1]] = list(mu=7)
init.cond1[[2]] = list(mu=7.5)
init.cond1[[3]] = list(mu=8)
init = init.cond1

FitNorm = "
model{
  for(i in 1:2){beta[i]~dnorm(0,.001)}       ## prior mean
  sigma ~ dgamma(0.001,0.001)   ## prior residual precision
  tau.t ~ dgamma(0.001,0.001)   ## prior year-effect precision
  
  for(t in 1:nt){           ## loop over years
    alpha.t[t] ~ dnorm(0,tau.t)     ## random year effect
    Ex[t] <- beta[1] + beta[2]*precip[t] + alpha.t[t]        ## process model 
    for(i in 1:nrep){           ## loop over reps
        
        x[t,i] ~ dnorm(Ex[t],sigma) ## data model
    }
    px[t] ~ dnorm(Ex[t],sigma)  ## predictive interval
  }
}
"
j.model   <- jags.model (file = textConnection(FitNorm),
                         data = data,
                         inits = NULL,
                         n.chains = 3)
b1   <- coda.samples (model = j.model,
                      variable.names = c("Ex","px","sigma","tau.t","alpha.t"),
                      n.iter = 130000,
                      burnin=2000,
                      thin=25)
bmcmc1 = b1

bmcmc1.df <- as.data.frame(as.matrix(bmcmc1))

quants <- rbind(
  quantile(bmcmc1.df$'Ex[1]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[2]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[3]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[4]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[5]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[6]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[7]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[8]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[9]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'Ex[10]',c(0.025,0.5,0.975)))

preds <- rbind(
  quantile(bmcmc1.df$'px[1]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[2]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[3]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[4]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[5]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[6]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[7]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[8]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[9]',c(0.025,0.5,0.975)),
  quantile(bmcmc1.df$'px[10]',c(0.025,0.5,0.975)))


ymin=min(preds[,1]);ymax=max(preds[,3])
plot(mosq$time,mosq$rep1, type="l", col=8, ylim=c(ymin, ymax))
lines(mosq$time,mosq$rep2, col=8)
lines(mosq$time,mosq$rep3, col=8)
lines(mosq$time,mosq$rep4, col=8)
lines(mosq$time,mosq$rep5, col=8)
lines(mosq$time,quants[,1],lty=2, col=2, lwd=2)
lines(mosq$time,quants[,2],col=2, lwd=2)
lines(mosq$time,quants[,3],lty=2,col=2, lwd=2)
lines(mosq$time,preds[,1],lty=2, col=3, lwd=2)
lines(mosq$time,preds[,3],lty=2, col=3, lwd=2)
