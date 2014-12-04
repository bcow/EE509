#!/usr/bin/env Rscript

library(MASS)

load("/usr2/collab/ecowdery/EE509/Project/mus.Rdata")

v1 = c(); v2 = c();
for(i in 1:5){
  v1 = c(v1, rep(i, 6-i))
  v2 = c(v2, (i+1):6)
}

k = nrow(out1.df)

tm1 <- system.time(
{
  DEN.1 <- list()
  for(N in 1:length(v1)){
    x <- out1.df[,v1[N]][1:k]
    y <- out1.df[,v2[N]][1:k]
    den <- kde2d(x, y, n = length(x))
    DEN.1 <- append(DEN.1, list(den))
    print(N)
  }
})
print(tm1)

save(DEN.1,file="/usr2/collab/ecowdery/EE509/Project/DEN1.Rdata")