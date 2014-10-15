

## to simulate data
x <- data.frame(onset=sample(as.Date("2014-01-01")+0:10, 30, replace=TRUE), patch=sample(c('a','b','c'), replace=TRUE, 30))




model <- function(x, w, D.patches=NULL, spa.kernel=dexp,
                  n.iter=1000, sample.every=200,
                  delta.ini=1,
                  phi.ini=0.001,
                  move.R=TRUE, sd.R=0.005,
                  move.delta=TRUE, sd.delta=0.001,
                  move.phi=TRUE, sd.phi=0.0001
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
    if(any(D.patches)<0) warning("D.patches has negative entries")
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
    colnames(incid) <- patches
    rownames(incid) <- as.character(all.dates)

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


    ## FUNCTIONS TO GET RATES AND LOG-LIKELIHOOD ##
    ## POISSON LIKELIHOOD
    ## p(incidence | rates)
    LL.all <- function(R, delta, phi){
        rates <- as.vector(get.lambdas(R, delta, phi))
        return(sum(dpois(incid.vec, rates, log=TRUE)))
    }


    ## PARAMETER MOVEMENTS ##
    ## MOVE R
    R.ACC <- 0
    R.REJ <- 0
    R.move <- function(R, sigma=sd.R){
        ## generate proposals ##
        newR <- R + rnorm(n=length(R), mean=0, sd=sd.R)

        if(all(newR>=0)){
            if((r <- log(runif(1))) <=  (LL.all(newR, delta, phi) - LL.all(R, delta, phi))){
                R <- temp # accept
                R.ACC <<- R.ACC+1
            } else { # reject
                R.REJ <<- R.REJ+1
            }
        } else { # reject
            R.REJ <<- R.REJ+1
        }

        ## return moved vector
        return(R)
    }


    ## MOVE SPATIAL PARAM 'DELTA'
    delta.ACC <- 0
    delta.REJ <- 0
    delta.move <- function(delta, sigma=sd.delta){
        ## generate proposals ##
        newdelta <- delta + rnorm(n=length(delta), mean=0, sd=sd.delta)

        if(all(newdelta>=0)){
            if((r <- log(runif(1))) <=  (LL.all(newdelta, delta, phi) - LL.all(delta, delta, phi))){
                delta <- temp # accept
                delta.ACC <<- delta.ACC+1
            } else { # reject
                delta.REJ <<- delta.REJ+1
            }
        } else { # reject
            delta.REJ <<- delta.REJ+1
        }

        ## return moved vector
        return(delta)
    }


    ## MOVE BACKGROUND FORCE OF INFECTION 'phi'
    phi.ACC <- 0
    phi.REJ <- 0
    phi.move <- function(phi, sigma=sd.phi){
        ## generate proposals ##
        newphi <- phi + rnorm(n=length(phi), mean=0, sd=sd.phi)

        if(all(newphi>=0)){
            if((r <- log(runif(1))) <=  (LL.all(newphi, phi, phi) - LL.all(phi, phi, phi))){
                phi <- temp # accept
                phi.ACC <<- phi.ACC+1
            } else { # reject
                phi.REJ <<- phi.REJ+1
            }
        } else { # reject
            phi.REJ <<- phi.REJ+1
        }

        ## return moved vector
        return(phi)
    }




    ## MCMC ##
    ## BASIC HEADER
    header <- "step\tpost\tlikelihood\tprior"

      cat(header, file=file.out)


    ## add first line
    ## temp: c(loglike, logprior)
    temp <- c(LL.all(y, x, phi, alpha), LPrior.alpha(alpha))

    ## check that initial LL is not -Inf
    if(!is.finite(temp[1])) warning("Initial likelihood is zero")

    ## write to file
    cat("\n", file=file.out, append=TRUE)
    cat(c(1, sum(temp), temp), sep="\t", append=TRUE, file=file.out)
    if(move.alpha) cat("", alpha/sum(alpha), sep="\t", append=TRUE, file=file.out)
    if(move.phi){
        cat("", as.vector(phi), sep="\t", append=TRUE, file=file.out)
    } else if(model.unsampled){
        cat("", as.vector(phi[,K]), sep="\t", append=TRUE, file=file.out)
    }

    if(!quiet) cat("\nStarting MCMC: 1")

    ## mcmc ##
    for(i in 1:n){
        ## move stuff ##
        ## move alpha if needed
        if(move.alpha) alpha <- alpha.move(alpha, sd.alpha)

        ## move phi if needed (phi.move makes the necessary moves)
        if((move.phi|model.unsampled) && (i %% move.phi.every == 0)) phi <- phi.move(phi, sd.phi)

        ## if retain this sample ##
        if(i %% sample.every ==0){
            temp <- c(LL.all(y, x, phi, alpha), LPrior.alpha(alpha))
            cat("\n", file=file.out, append=TRUE)
            cat(c(i, sum(temp), temp), sep="\t", append=TRUE, file=file.out)
            if(move.alpha) cat("", alpha/sum(alpha), sep="\t", append=TRUE, file=file.out)
            if(move.phi) {
                cat("", as.vector(phi), sep="\t", append=TRUE, file=file.out)
            } else if(model.unsampled){
                cat("", as.vector(phi[,K]), sep="\t", append=TRUE, file=file.out)
            }
            if(!quiet) cat("..",i)
        }
    }

    if(!quiet) cat("..done!\nResults were saved in file:",file.out,"\n")


    ## re-read output file ##
    out <- read.table(file.out, header=TRUE, colClasses="numeric", sep="\t")
    out$step <- as.integer(out$step)

    ## format using coda ##
    if(!quiet){
        ## acceptance rates for alpha
        if(move.alpha){
            cat("\nacceptance rate for alpha: ", ALPHA.ACC/(ALPHA.ACC+ALPHA.REJ))
            cat("\naccepted: ", ALPHA.ACC)
            cat("\nreject: ", ALPHA.REJ)
            cat("\n")
        }

        ## acceptance rates for phi
        if(move.phi || model.unsampled){
            cat("\nacceptance rate for phi: ", PHI.ACC/(PHI.ACC+PHI.REJ))
            cat("\naccepted: ", PHI.ACC)
            cat("\nreject: ", PHI.REJ)
            cat("\n")
        }

    }


    class(out) <- c("data.frame", "bmmix")
    return(out)

 
    

}
