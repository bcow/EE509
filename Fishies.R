x = c(5,6,7,3,6,5,8,4,4,3)
hist(x,c(.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5))
lines(x,dpois(x,1.96),type='s')
lines(dpois(x,2),col=2)

















lines(density(x))

d <- density(x)
plot(d)


lklh.poisson <- function(x, lambda) lambda^x/factorial(x) * exp(-lambda)

log.lklh.poisson <- function(x, lambda){ 
    -sum(x * log(lambda) - log(factorial(x)) - lambda) 
  }


optim(par = 2, log.lklh.poisson, x = x)

fitdistr(x,densfun="negative binomial")
sum(x)
10/51

y=seq(0,9,.01)
lines(y,dpois(y,2),type='s')

plot(x,dpois(x,2),type='s')
