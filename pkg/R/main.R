

## to simulate data
x <- data.frame(onset=sample(as.Date("2014-01-01")+0:10, 30, replace=TRUE), patch=sample(c('a','b','c'), replace=TRUE, 30))




epidemicMCMC <- function(x, w, D.patches=NULL, spa.kernel=dexp,
                         n.iter=1e5, sample.every=200,
                         move.R=TRUE, sd.R=0.005, R.ini=1,
                         move.delta=TRUE, sd.delta=0.001, delta.ini=1,
                         move.phi=TRUE, sd.phi=0.0001, phi.ini=0.001,
                         tune=TRUE, max.tune=2e4,
                         file.out="mcmc.txt", quiet=FALSE){

    ## CHECKS / HANDLE ARGUMENTS ##
    x <- na.omit(x)
    x$patch <- factor(x$patch)
    patches <- levels(x$patch)
    n.patches <- length(patches)

    ## handle empty distance matrix between patches
    if(is.null(D.patches)){
        D.patches <- matrix(0, ncol=n.patches,nrow=n.patches)
        colnames(D.patches) <- rownames(D.patches) <- patches
    }

    ## force null diagonal
    diag(D.patches) <- 0

    ## useful warnings/errors
    if(any(D.patches<0)) warning("D.patches has negative entries")
    if(ncol(D.patches)!= nrow(D.patches)) stop(paste("D.patches is not square",
           nrow(D.patches),"rows,",
           ncol(D.patches), "columns."))
    if(ncol(D.patches)!= n.patches) stop(paste("D.patches has wrong dimensions: expected",
           n.patches, "but got",
           ncol(D.patches)))


    ## HANDLE DATES
    date.range <- range(x$onset)
    all.dates <- seq(date.range[1],date.range[2],by=1)
    n.dates <- length(all.dates)

    ## HANDLE SERIAL INTERVAL DISTRIBUTION
    ## note: w[i] is the proba of infection after i days
    w <- w/sum(w) # standardize distri to 1
    w <- c(w, rep(0, n.dates))


    ## COMPUTE INCIDENCE ##
    ## function to get daily incidence
    get.daily.incid <- function(dates, from, to){
        out <- table(factor(as.character(dates), levels=as.character(seq(from, to, by=1))))
        return(out)
    }

    ## daily incidence
    temp <- tapply(x$onset, x$patch, get.daily.incid,
                   from=date.range[1],
                   to=date.range[2])

    ## put incidence in the right format
    ## basic matrix
    incid.mat <- Reduce(cbind,temp)
    colnames(incid.mat) <- patches
    rownames(incid.mat) <- as.character(all.dates)

    ## store as vector for faster computations in 'dpois'
    incid.vec <- as.vector(incid.mat)


    ## COMPUTE INFECTIOUSNESS (CONSTANT PART OF IT) ##
    ## AUXILIARY FUNCTION:
    ## gets infectiousness for day 'date' and
    ## a vector of onset dates 'onsets'
    get.infec <- function(date, onsets){
        temp <- as.integer(date-onsets)
        out <- sum(w[temp[temp>0]])
        return(out)
    }

    ## GET INFECTIOUSNESS
    ## 'infec.mat' contains the sum of infectiousness for each day (row) and patch (column)
    ## basic matrix
    infec.mat <- t(sapply(all.dates, function(t) tapply(x$onset, x$patch, function(e) get.infec(t, e))))
    infec.per.patch <- apply(infec.mat,2,sum)

    ## GET BETA FOR ALL PATCHES ##
    ## (force of infection coming from a given patch at a given time)
    ## note: R can be a vector
    get.betas <- function(R){
        return(R*infec.per.patch)
    }


    ## GET LAMBDA FOR ALL PATCHES ##
    ## (force of infection experienced by the patches)
    ## computed using matrix product
    ## (beta %*% P) + phi, with:
    ## - beta: infectiousness for day (row) x patches (col)
    ## - P: connectivity from patches (row) to patches (col),
    ## standardized by column
    ## - phi: the background force of infection
    ## here, beta = R * infect.mat
    ## and P is given by 'spa.kernel(D.pacthes, delta)'
    ## but needs to be column-standardized
    ##
    get.lambdas <- function(R, delta, phi){
        return(
            (R*infec.mat %*%
             prop.table(spa.kernel(D.patches, delta),2)
             ) + phi
            )
    }


    ## LOG-LIKELIHOOD AND PRIORS ##
    ## POISSON LIKELIHOOD
    ## p(incidence | rates)
    LL <- function(R, delta, phi){
        rates <- as.vector(get.lambdas(R, delta, phi))
        return(sum(dpois(incid.vec, rates, log=TRUE)))
    }

    ## PRIORS ##
    ## prior for R
    R.logprior <- function(R) return(0)

    ## prior for delta
    delta.logprior <- function(delta) return(0)

    ## prior for phi
    phi.logprior <- function(phi) return(0)

    ## all priors
    all.logprior <- function(R, delta, phi){
        return(R.logprior(R) + delta.logprior(delta) + phi.logprior(phi))
    }



    ## PARAMETER MOVEMENTS ##
    ## MOVE R
    R.ACC <- 0
    R.REJ <- 0
    R.move <- function(R, delta, phi, sigma=sd.R){
        ## generate proposals ##
        newR <- R + rnorm(n=length(R), mean=0, sd=sd.R)

        if(all(newR>=0)){
            if((r <- log(runif(1))) <=  (LL(newR, delta, phi) - LL(R, delta, phi) +
                                         R.logprior(newR) - R.logprior(R))){
                R <- newR # accept
                R.ACC <<- R.ACC+1
            } else { # reject
                R.REJ <<- R.REJ+1
            }
        } else { # reject
            R.REJ <<- R.REJ+1
        }

        ## return moved vector
        return(R)
    } # end R.move


    ## MOVE SPATIAL PARAM 'DELTA'
    delta.ACC <- 0
    delta.REJ <- 0
    delta.move <- function(R, delta, phi, sigma=sd.delta){
        ## generate proposals ##
        newdelta <- delta + rnorm(n=length(delta), mean=0, sd=sd.delta)

        if(all(newdelta>=0)){
            if((r <- log(runif(1))) <=  (LL(R, newdelta, phi) - LL(R, delta, phi) +
                                         delta.logprior(newdelta) - delta.logprior(delta))){
                delta <- newdelta # accept
                delta.ACC <<- delta.ACC+1
            } else { # reject
                delta.REJ <<- delta.REJ+1
            }
        } else { # reject
            delta.REJ <<- delta.REJ+1
        }

        ## return moved vector
        return(delta)
    } # end delta.move


    ## MOVE BACKGROUND FORCE OF INFECTION 'phi'
    phi.ACC <- 0
    phi.REJ <- 0
    phi.move <- function(R, delta, phi, sigma=sd.phi){
        ## generate proposals ##
        newphi <- phi + rnorm(n=length(phi), mean=0, sd=sd.phi)

        if(all(newphi>=0)){
            if((r <- log(runif(1))) <=  (LL(R, delta, newphi) - LL(R, delta, phi) +
                                         phi.logprior(newphi) - phi.logprior(phi))){
                phi <- newphi # accept
                phi.ACC <<- phi.ACC+1
            } else { # reject
                phi.REJ <<- phi.REJ+1
            }
        } else { # reject
            phi.REJ <<- phi.REJ+1
        }

        ## return moved vector
        return(phi)
    } # end phi.move




    ## TUNING ##
    stop.tune.R <- stop.tune.delta <- stop.tune.phi <- FALSE

    ## TUNING FUNCTION FOR R
    R.tune <- function(sd.R){
        ## jumps too large
        if((R.ACC/(R.ACC + R.REJ))<0.2){
            return(sd.R*0.8)
        }

        ## jumps too small
        if((R.ACC/(R.ACC + R.REJ))>0.5){
            return(sd.R*1.2)
        }

        ## if reaching here, stop tuning R
        stop.tune.R <<- TRUE
    } # end R.tune


    ## TUNING FUNCTION FOR DELTA
    delta.tune <- function(sd.delta){
        ## jumps too large
        if((delta.ACC/(delta.ACC + delta.REJ))<0.2){
            return(sd.delta*0.8)
        }

        ## jumps too small
        if((delta.ACC/(delta.ACC + delta.REJ))>0.5){
            return(sd.delta*1.2)
        }

        ## if reaching here, stop tuning delta
        stop.tune.delta <<- TRUE
    } # end delta.tune


    ## TUNING FUNCTION FOR PHI
    phi.tune <- function(sd.phi){
        ## jumps too large
        if((phi.ACC/(phi.ACC + phi.REJ))<0.2){
            return(sd.phi*0.8)
        }

        ## jumps too small
        if((phi.ACC/(phi.ACC + phi.REJ))>0.5){
            return(sd.phi*1.2)
        }

        ## if reaching here, stop tuning phi
        stop.tune.phi <<- TRUE
    } # end phi.tune


    ## INITIALIZE VALUES
    R <- R.ini
    delta <- delta.ini
    phi <- phi.ini

    COUNTER <- 0
    while(tune && (COUNTER<=max.tune)){
        ## update counter
        COUNTER <- COUNTER + 1

        ## move stuff ##
        ## move R if needed
        if(move.R) R <- R.move(R, delta, phi)

        ## move delta if needed
        if(move.delta) delta <- delta.move(R, delta, phi)

        ## move phi if needed
        if(move.phi) phi <- phi.move(R, delta, phi)

        ## change parameters of proposal distributions ##
        ## (every 100 iterations)
        if(COUNTER %% 100 ==0){
            ## update SDs
            sd.R <- R.tune(sd.R)
            sd.delta <- delta.tune(sd.delta)
            sd.phi <- phi.tune(sd.phi)

            ## check if we should stop
            tune <- !(all(stop.tune.R, stop.tune.delta, stop.tune.phi))
        }
    }



    ## MCMC ##
    ## INITIALIZE VALUES
    R <- R.ini
    delta <- delta.ini
    phi <- phi.ini

    ## BASIC HEADER
    header <- "step\tpost\tlikelihood\tprior\tR\tdelta\tphi"
    cat(header, file=file.out)

    ## add first line
    ## temp: c(loglike, logprior)
    temp <- c(LL(R, delta, phi), all.logprior(R, delta, phi))

    ## check that initial LL is not -Inf
    if(!is.finite(temp[1])) warning("Initial likelihood is zero")

    ## write to file
    cat("\n", file=file.out, append=TRUE)
    cat(c(1, sum(temp), temp, R, delta, phi), sep="\t", append=TRUE, file=file.out)

    ## basic message
    if(!quiet) cat("\nStarting MCMC: 1")

    ## mcmc ##
    for(i in 1:n.iter){
        ## move stuff ##
        ## move R if needed
        if(move.R) R <- R.move(R, delta, phi)

        ## move delta if needed
        if(move.delta) delta <- delta.move(R, delta, phi)

        ## move phi if needed
        if(move.phi) phi <- phi.move(R, delta, phi)

        ## if retain this sample ##
        if(i %% sample.every ==0){
            temp <- c(LL(R, delta, phi), all.logprior(R, delta, phi))
            cat("\n", file=file.out, append=TRUE)

            ## add posterior
            cat(c(i, sum(temp), temp, R, delta, phi), sep="\t", append=TRUE, file=file.out)

            ## basic message
            if(!quiet) cat("..",i)
        }
    }

    ## basic message
    if(!quiet) cat("..done!\nResults were saved in file:",file.out,"\n")

    ## re-read output file ##
    out <- read.table(file.out, header=TRUE, colClasses="numeric", sep="\t")
    out$step <- as.integer(out$step)

    ## format using coda ##
    if(!quiet){
        ## acceptance rates for R
        if(move.R){
            cat("\nacceptance rate for R: ", R.ACC/(R.ACC+R.REJ))
            cat("\naccepted: ", R.ACC)
            cat("\nreject: ", R.REJ)
            cat("\n")
        }

        ## acceptance rates for delta
        if(move.delta){
            cat("\nacceptance rate for delta: ", delta.ACC/(delta.ACC+delta.REJ))
            cat("\naccepted: ", delta.ACC)
            cat("\nreject: ", delta.REJ)
            cat("\n")
        }

        ## acceptance rates for phi
        if(move.phi){
            cat("\nacceptance rate for phi: ", phi.ACC/(phi.ACC+phi.REJ))
            cat("\naccepted: ", phi.ACC)
            cat("\nreject: ", phi.REJ)
            cat("\n")
        }
    }


    class(out) <- c("data.frame", "epidemicMCMC")
    return(out)

} # end epidemicMCMC
