
## clean environment
rm(list=ls())

## source code
library(ggplot2)
source("../pkg/R/main.R")
source("../pkg/R/plotmcmc.R")
## source("~/dev/spatialPoisson/code/pkg/R/main.R")
## source("~/dev/spatialPoisson/code/pkg/R/plotmcmc.R")



## LOAD DATA
dat <- read.csv("WestAfrica.infevents.csv")

## PUT DATA IN RIGHT SHAPE
onset <- dat$t
patch <- dat$adunit_infectee
meanX <- tapply(dat$x_infectee, dat$adunit_infectee, mean)
meanY <- tapply(dat$y_infectee, dat$adunit_infectee, mean)
patch.xy <- cbind(meanX, meanY)
D.patches <- dist(patch.xy)

## check data
f.dat <- ggplot(data.frame(onset,patch=factor(patch)),aes(x=onset)) +
    geom_histogram(aes(fill=patch),bin=1) +
    labs(x="Date of onset",y="Number of new daily cases", title="Incidence")
f.dat


## get generation time
SI.gammaPar <- .gammaToMeanVar(14.2, 9.5)
SI <- dgamma(1:10, shape=SI.gammaPar[1], rate=SI.gammaPar[1])
barplot(SI,names=1:10, xlab="Days after onset", main="Serial interval")


res <- epidemicMCMC(onset[onset>20], patch[onset>20], w=SI)

par(mfrow=c(3,3))
for(i in 2:9) plot(res$chains[-1,1], res$chains[-1,i], main=names(res$chains)[i],type="l",xlab="",ylab="")


#### OLD DUMMY DATASETS ####
## simulate toy data
## set.seed(1)
## x <- data.frame(onset=sample(as.Date("2014-01-01")+0:10, 200, replace=TRUE), patch=sample(c('a','b','c','d','e'), replace=TRUE, 200))


## run MCMC
#set.seed(1)
#system.time(res <- epidemicMCMC(x,w=c(1,2,1),n.iter=5e4,sample.every=200, max.tune=1e4))

## system.time(res <- epidemicMCMC(x,w=c(1,2,1),n.iter=2e6,sample.every=500, max.tune=2e4))

## ## make some graphs
## pdf("alltraces.pdf")

## par(mfrow=c(3,3))
## for(i in 2:9) plot(res$chains[-1,1], res$chains[-1,i], main=names(res$chains)[i],type="l",xlab="",ylab="")

## dev.off()


## pdf("toyexpl-logpost.pdf")
## plotmcmc(res, what=c("post"))$traces
## dev.off()

## pdf("toyexpl-R.pdf")
## plotmcmc(res, what=c("R"))$traces
## dev.off()

## pdf("toyexpl-R.pdf")
## plotmcmc(res, what=c("R"))$densities
## dev.off()


