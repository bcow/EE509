alcesdata <- read.table("~/R Projects/EE509/data/alcesdata.txt", header=TRUE, quote="\"")

FitNorm = "
model{
  r~dnorm(0,.001)
  tau~dgamma(.001,.001)
  sigma~dgamma(.001,.001)

  x[1]~dnorm(0,.001) 
  y[1]~dnorm(x[1],tau)
  
  for (i in 2:N){ 
    u[i]<-x[i-1] + r
    x[i]~dnorm(u[i],sigma)
    y[i]~dnorm(x[i],tau)
  }
}
"
data = list(y=log(alcesdata$density), N=length(alcesdata$density))
init = NULL

j.model   <- jags.model (file = textConnection(FitNorm),
                         data = data,
                         inits = init,
                         n.chains = 3)
b1   <- coda.samples (model = j.model,
                      variable.names = c("r","tau","sigma",'y'),
                      n.iter=300000,
                      burnin=2000,
                      thin=50
                      )
bmcmc1 = b1
gelman.plot(bmcmc1) 
# autocorr.plot(bmcmc1)
cumuplot(bmcmc1)
plot(bmcmc1)
summary(bmcmc1)


# 1.  Write out and run the BUGS code for the exponential growth model. 
# a) Plots of the time series data on both the LOG and LINEAR scales that include the model mean and credible intervals.

par(mfrow=c(2,1))
plot(alcesdata$year,alcesdata$density)
plot(alcesdata$year,alcesdata$density,log='y')

# b) Density plots for the predicted values for the missing data

idx = which(is.na(alcesdata$density))
bmcmc1.df <- as.data.frame(as.matrix(bmcmc1))
idx

par(mfrow=c(3,2))
plot(density(bmcmc1.df$'y[18]'))
plot(density(bmcmc1.df$'y[19]'))
plot(density(bmcmc1.df$'y[20]'))
plot(density(bmcmc1.df$'y[21]'))
plot(density(bmcmc1.df$'y[44]'))

# c) Density plots of the intrinsic growth rate, the observation error variance, and the process model error variance

par(mfrow=c(1,3))
plot(density(bmcmc1.df$r))
plot(density(bmcmc1.df$sigma))
plot(density(bmcmc1.df$tau))

# 2.  Modify the exponential growth process model in the BUGS code to instead be the Ricker growth model.  Rerun including your BUGS code and the same figures as in part 1 plus plots for both the prior and posterior density on the carrying capacity.

FitNorm = "
model{
  r~dnorm(0,.001)
  K~dlnorm(0,.001)
  tau~dgamma(.001,.001)
  sigma~dgamma(.001,.001)

  x[1]~dnorm(0,.001) 
  y[1]~dnorm(x[1],tau)
  
  for (i in 2:N){ 
    u[i]<-x[i-1] + r*(1-exp(x[i-1])/K)

    x[i]~dnorm(u[i],sigma)
    y[i]~dnorm(x[i],tau)
  }
}
"
data = list(y=log(alcesdata$density), N=length(alcesdata$density))
init.cond1 <- list()
init.cond1[[1]] = list(k=250)
init.cond1[[2]] = list(k=255)
init.cond1[[3]] = list(k=260)
init = NULL

j.model   <- jags.model (file = textConnection(FitNorm),
                         data = data,
                         inits = init,
                         n.chains = 3)
update(j.model, 2000)
b1   <- coda.samples (model = j.model,
                      variable.names = c("r",'K',"tau","sigma"),
                      n.iter=100000,
                      thin=50
)
bmcmc1 = b1
gelman.plot(bmcmc1) 
autocorr.plot(bmcmc1)
cumuplot(bmcmc1)
plot(bmcmc1)
summary(bmcmc1)


plot(density(bmcmc1.df$K))


# r = beta[1]
# beta[2] = r/K --> K= r/beta[1]
# 
# beta[1]+beta[2]*exp(x[i])

FitNorm = "
model{

  for(i in 1:2){beta[i]~dnorm(0,.001)}
  tau~dgamma(.001,.001)
  sigma~dgamma(.001,.001)

  x[1]~dnorm(0,.001) 
  y[1]~dnorm(x[1],tau)

  for (i in 2:N){ 
    u[i]<-x[i-1] + beta[1] + beta[2]*exp(x[i-1])
    x[i]~dnorm(u[i],sigma)
    y[i]~dnorm(x[i],tau)

  }
}
"
j.model   <- jags.model (file = textConnection(FitNorm),
                         data = data,
                         inits = init,
                         n.chains = 3)
b2   <- coda.samples (model = j.model,
                      variable.names = c("beta","tau","sigma"),
                      n.iter=100000,
                      burnin=2000,
                      thin=50
                      )
bmcmc2 = b2
summary(bmcmc2)
