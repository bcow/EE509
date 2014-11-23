data <- list( nu=5, alpha=c(370,370,370),
      Tau = structure(.Data = c(0.000001, 0, 0, 0, 0.000001, 0, 0, 0, 0.000001), .Dim = c(3, 3)),
      Lambda = structure(.Data = c(2000, 0, 0, 0, 2000, 0, 0, 0, 2000), .Dim = c(3, 3)),
      Y = structure(.Data = c(
        343,463,NA,
        332,478,NA,
        365,378,NA,
        332,448,NA,
        353,465,NA,
        357,546,NA,
        346,473,NA,
        343,461,NA,
        378,NA,381,
        345,NA,372,
        354,NA,365,
        347,NA,362,
        351,NA,349,
        340,NA,330,
        365,NA,365,
        345,NA,338,
        360,NA,382,
        357,NA,383,
        345,NA,376,
        347,NA,370,
        NA,495,NA,
        NA,509,NA,
        NA,521,NA,
        NA,529,NA,
        NA,568,NA,
        NA,554,NA,
        NA,494,NA,
        NA,573,NA,
        NA,523,NA,
        NA,496,NA,
        NA,535,NA,
        NA,509,NA,
        NA,442,NA,
        NA,481,NA,
        NA,481,NA,
        NA,468,NA,
        NA,NA,333,
        NA,NA,336,
        NA,NA,331,
        NA,NA,320,
        NA,NA,361,
        NA,NA,370,
        NA,NA,387,
        NA,NA,385,
        NA,NA,342,
        NA,NA,364,
        NA,NA,315,
        NA,NA,430,
        NA,NA,365,
        NA,NA,338,
        NA,NA,382), .Dim = c(51, 3)))

Model = "
model{
  for(i in 1:51){
    Y[i, 1:3] ~ dmnorm(mu[], R[ , ])
  }
  mu[1:3] ~ dmnorm(alpha[],Tau[ , ])
  
  R[1:3 , 1:3] ~ dwish(Lambda[ , ], nu)
  D[1:3, 1:3]<-inverse(R[1:3, 1:3])
  sig1<-sqrt(D[1,1])
  sig2<-sqrt(D[2,2])
  sig3<-sqrt(D[3,3])
  
  rho12<-D[1,2]/(sig1*sig2)
  rho13<-D[1,3]/(sig1*sig3)
  rho23<-D[2,3]/(sig2*sig3)
  
  diff21<-mu[2]-mu[1]
  diff31<-mu[3]-mu[1]
}"

init = NULL

j.model   <- jags.model (file = textConnection(Model),
                         data = data,
                         inits = init,
                         n.chains = 3)