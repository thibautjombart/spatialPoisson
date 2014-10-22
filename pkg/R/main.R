

## to simulate data
x <- data.frame(onset=sample(as.Date("2014-01-01")+0:300, 5000, replace=TRUE), patch=sample(c('a','b','c','d','e'), replace=TRUE, 5000))



##################
## epidemicMCMC ##
##################
## Main function implementing the patch Poisson model
## Arguments
## - x: a data.frame with column 'onset' (Date object), and column 'patch'
## - w: a numeric describing the generation time / serial interval distribution
## (w[1] is the proba of infection 1 day after infection / symptoms)
## - D.patches: square matrix of distances between patches
## - spa.kernel: a function describing the spatial kernel
## - n.iter: length of the MCMC
## - sample.every: frequency at which chains are saved
## - move.[param]: a logical indicating if 'param' should be moved
## - sd.[param]: the standard deviation of the proposal distribution for 'param'
## - [param].ini: initial value of 'param'
## - prior.[param]: parameters of the prior for 'param'
## - tune: a logical indicating if auto-tuning should be done (target: 20-50% acceptance)
## - max.tune: the maximum number of iterations for the auto-tuning
## - file.out: path to the file where MCMC output are written
## - quiet: a logical indicating if messages should be hidden
epidemicMCMC <- function(x, w, D.patches=NULL, spa.kernel=dexp,
                         n.iter=1e5, sample.every=200,
                         move.R=TRUE, sd.R=0.005, R.ini=1,
                         move.delta=TRUE, sd.delta=0.001, delta.ini=1,
                         move.rho=TRUE, sd.rho=0.0001, rho.ini=0.001,
                         logprior.delta=function(x) dexp(x, rate=1,log=TRUE),
                         logprior.rho=function(x) 0,
                         tune=TRUE, max.tune=2e4,
                         file.out="mcmc.txt", quiet=FALSE){

    ## CHECKS / HANDLE ARGUMENTS ##
    x <- na.omit(x)
    x$patch <- factor(x$patch)
    patches <- levels(x$patch)
    n.patches <- length(patches)
    x$day <- as.integer(x$onset-min(x$onset))

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


    ## GET INFECTIOUSNESS ##

    ## GET BETA FOR ALL PATCHES AT TIME T
    ## beta_t_i: infectiousness due to individuals in patch 'i' at time 't'
    ## beta_t^i = sum_{s=1}^{t-1} N_s^i w(t-s)
    ## N is a matrix of incidence (row=dates,col=patches)
    get.betas.t <- function(N,t) return(apply(incid.mat[1:(t-1),,drop=FALSE] * w[(t-1):1], 2, sum))

    ## GET ALL BETAS (excluding t=1)
    get.betas <- function(N) {
        return(t(sapply(2:nrow(N), function(t) get.beta.t(N,t))))
    }

    ## ## GET ALL BETAS (second version, scales poorly with number of cases)
    ## get.beta.all <- function(x) tapply(x$day, x$patch, function(dates) # do by patch
    ##                                    sapply(2:nrow(x), function(t)sum(w[t-dates[dates<t]]))) # sum infectiousnesses



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
    ## (beta %*% P) with:
    ## - beta: infectiousness for day (row) x patches (col)
    ## - P: connectivity from patches (row) to patches (col),
    ## standardized by column
    ## - rho: the background force of infection
    ## here, beta = R * infect.mat
    ## and P is given by 'spa.kernel(D.pacthes, delta)'
    ## but needs to be column-standardized
    ##
    get.lambdas <- function(R, delta, rho){
        return(
            (R*infec.mat %*%
             prop.table(spa.kernel(D.patches, delta),2)
             )
            )
    }


    ## LOG-LIKELIHOOD AND PRIORS ##
    ## POISSON LIKELIHOOD
    ## p(incidence | rates)
    LL <- function(R, delta, rho){
        rates <- as.vector(get.lambdas(R, delta, rho))
        return(sum(dpois(incid.vec, rates, log=TRUE)))
    }

    ## p(R|rho)
    LL.R <- function(R,rho){
        return(dgamma(R, shape=rho[1], rate=rho[2], log=TRUE))
    }

    ## PRIORS ##
    ## all priors
    logprior.all <- function(delta, rho){
        return(logprior.delta(delta) + logprior.rho(rho))
    }



    ## PARAMETER MOVEMENTS ##
    ## MOVE R
    R.ACC <- 0
    R.REJ <- 0
    R.move <- function(R, delta, rho, sigma=sd.R){
        ## generate proposals ##
        newR <- R + rnorm(n=length(R), mean=0, sd=sd.R)

        if(all(newR>=0)){
            if((r <- log(runif(1))) <=  (LL(newR, delta, rho) - LL(R, delta, rho) +
                                         logprior.R(newR) - logprior.R(R))){
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
    delta.move <- function(R, delta, rho, sigma=sd.delta){
        ## generate proposals ##
        newdelta <- delta + rnorm(n=length(delta), mean=0, sd=sd.delta)

        if(all(newdelta>=0)){
            if((r <- log(runif(1))) <=  (LL(R, newdelta, rho) - LL(R, delta, rho) +
                                         logprior.delta(newdelta) - logprior.delta(delta))){
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


    ## MOVE BACKGROUND FORCE OF INFECTION 'rho'
    rho.ACC <- 0
    rho.REJ <- 0
    rho.move <- function(R, delta, rho, sigma=sd.rho){
        ## generate proposals ##
        newrho <- rho + rnorm(n=length(rho), mean=0, sd=sd.rho)

        if(all(newrho>=0)){
            if((r <- log(runif(1))) <=  (LL(R, delta, newrho) - LL(R, delta, rho) +
                                         logprior.rho(newrho) - logprior.rho(rho))){
                rho <- newrho # accept
                rho.ACC <<- rho.ACC+1
            } else { # reject
                rho.REJ <<- rho.REJ+1
            }
        } else { # reject
            rho.REJ <<- rho.REJ+1
        }

        ## return moved vector
        return(rho)
    } # end rho.move




    ## TUNING ##
    stop.tune.R <- stop.tune.delta <- stop.tune.rho <- FALSE

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


    ## TUNING FUNCTION FOR RHO
    rho.tune <- function(sd.rho){
        ## jumps too large
        if((rho.ACC/(rho.ACC + rho.REJ))<0.2){
            return(sd.rho*0.8)
        }

        ## jumps too small
        if((rho.ACC/(rho.ACC + rho.REJ))>0.5){
            return(sd.rho*1.2)
        }

        ## if reaching here, stop tuning rho
        stop.tune.rho <<- TRUE
    } # end rho.tune


    ## INITIALIZE VALUES
    R <- R.ini
    delta <- delta.ini
    rho <- rho.ini

    ## basic message
    if(!quiet && tune) cat("Starting tuning proposal distributions...\n")

    COUNTER <- 0
    KEEPTUNING <- tune
    while(KEEPTUNING && (COUNTER<=max.tune)){
        ## update counter
        COUNTER <- COUNTER + 1

        ## move stuff ##
        ## move R if needed
        if(move.R) R <- R.move(R, delta, rho)

        ## move delta if needed
        if(move.delta) delta <- delta.move(R, delta, rho)

        ## move rho if needed
        if(move.rho) rho <- rho.move(R, delta, rho)

        ## change parameters of proposal distributions ##
        ## (every 100 iterations)
        if(COUNTER %% 100 ==0){
            ## update SDs
            if(move.R) sd.R <- R.tune(sd.R) else stop.tune.R <- TRUE
            if(move.delta) sd.delta <- delta.tune(sd.delta) else stop.tune.delta <- TRUE
            if(move.rho) sd.rho <- rho.tune(sd.rho) else stop.tune.rho <- TRUE

            ## check if we should stop
            KEEPTUNING <- !(all(stop.tune.R, stop.tune.delta, stop.tune.rho))
        }
    }

    ## basic message
    if(!quiet && tune){
        cat("... tuning done\n")
        ## acceptance rates for R
        if(move.R){
            cat("\ntuned acceptance rate for R: ", R.ACC/(R.ACC+R.REJ))
            cat("\naccepted: ", R.ACC)
            cat("\nreject: ", R.REJ)
            cat("\ntotal: ", R.ACC + R.REJ)
            cat("\n")
        }

        ## acceptance rates for delta
        if(move.delta){
            cat("\ntuned acceptance rate for delta: ", delta.ACC/(delta.ACC+delta.REJ))
            cat("\naccepted: ", delta.ACC)
            cat("\nreject: ", delta.REJ)
            cat("\ntotal: ", delta.ACC + delta.REJ)
            cat("\n")
        }

        ## acceptance rates for rho
        if(move.rho){
            cat("\ntuned acceptance rate for rho: ", rho.ACC/(rho.ACC+rho.REJ))
            cat("\naccepted: ", rho.ACC)
            cat("\nreject: ", rho.REJ)
            cat("\ntotal: ", rho.ACC + rho.REJ)
            cat("\n")
        }
    }



    ## MCMC ##
    ## INITIALIZE VALUES
    R <- R.ini
    delta <- delta.ini
    rho <- rho.ini

    ## BASIC HEADER
    header <- "step\tpost\tlikelihood\tprior\tR\tdelta\trho"
    cat(header, file=file.out)

    ## add first line
    ## temp: c(loglike, logprior)
    temp <- c(LL(R, delta, rho), logprior.all(R, delta, rho))

    ## check that initial LL is not -Inf
    if(!is.finite(temp[1])) warning("Initial likelihood is zero")

    ## write to file
    cat("\n", file=file.out, append=TRUE)
    cat(c(1, sum(temp), temp, R, delta, rho), sep="\t", append=TRUE, file=file.out)

    ## basic message
    if(!quiet) cat("\nStarting MCMC: 1")

    ## mcmc ##
    for(i in 1:n.iter){
        ## move stuff ##
        ## move R if needed
        if(move.R) R <- R.move(R, delta, rho)

        ## move delta if needed
        if(move.delta) delta <- delta.move(R, delta, rho)

        ## move rho if needed
        if(move.rho) rho <- rho.move(R, delta, rho)

        ## if retain this sample ##
        if(i %% sample.every ==0){
            temp <- c(LL(R, delta, rho), logprior.all(R, delta, rho))
            cat("\n", file=file.out, append=TRUE)

            ## add posterior
            cat(c(i, sum(temp), temp, R, delta, rho), sep="\t", append=TRUE, file=file.out)

            ## basic message
            if(!quiet) cat("..",i)
        }
    }

    ## basic message
    if(!quiet) cat("..done!\nResults were saved in file:",file.out,"\n")

    ## BUILD OUTPUT ##
    out <- list()
    ## re-read output file ##
    out$chains <- read.table(file.out, header=TRUE, colClasses="numeric", sep="\t")
    out$chains$step <- as.integer(out$chains$step)

    ## matched call ##
    out$call <- match.call()

    ## acceptance rates ##
    out$AR <- data.frame(accept=c(R.ACC,delta.ACC,rho.ACC),
               reject=c(R.REJ,delta.REJ,rho.REJ))
    out$AR$prop <- out$AR$accept / (out$AR$accept + out$AR$reject)
    rownames(out$AR) <- c("R","delta","rho")

    ## add class, and return ##
    class(out) <- c("list", "epidemicMCMC")
    return(out)

} # end epidemicMCMC





