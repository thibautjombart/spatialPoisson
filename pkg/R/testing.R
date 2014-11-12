
## source code
library(ggplot2)
source("~/dev/spatialPoisson/code/pkg/R/main.R")
source("~/dev/spatialPoisson/code/pkg/R/plotmcmc.R")

## simulate toy data
set.seed(1)
x <- data.frame(onset=sample(as.Date("2014-01-01")+0:100, 2000, replace=TRUE), patch=sample(c('a','b','c','d','e'), replace=TRUE, 2000))


## run MCMC
set.seed(1)
system.time(res <- epidemicMCMC(x,w=c(1,2,1),n.iter=5e4,sample.every=200, max.tune=1e4))

system.time(res <- epidemicMCMC(x,w=c(1,2,1),n.iter=1e5,sample.every=500))

## make some graphs
pdf("alltraces.pdf")

par(mfrow=c(3,3))
for(i in 2:9) plot(res$chains[-1,1], res$chains[-1,i], main=names(res$chains)[i],type="l",xlab="",ylab="")

dev.off()


## pdf("toyexpl-logpost.pdf")
## plotmcmc(res, what=c("post"))$traces
## dev.off()

## pdf("toyexpl-R.pdf")
## plotmcmc(res, what=c("R"))$traces
## dev.off()

## pdf("toyexpl-R.pdf")
## plotmcmc(res, what=c("R"))$densities
## dev.off()


