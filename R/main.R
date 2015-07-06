
## to simulate data
## x <- data.frame(onset=sample(as.Date("2014-01-01")+0:300, 5000, replace=TRUE), patch=sample(c('a','b','c','d','e'), replace=TRUE, 5000))
## w <- c(1,2,1)
## D.patches=NULL

spa.kernel=dexp
n.iter=1e5
sample.every=200
move.R=TRUE
sd.R=0.01
R.ini=1
move.delta=TRUE
sd.delta=0.001
delta.ini=1
move.pi=TRUE
sd.pi=0.001
pi.ini=0.5
move.N=TRUE
N.ini=NULL
sd.N=1
move.N.every=50
logprior.R=function(x) dgamma(R, shape=2, rate=1, log=TRUE)
logprior.delta=function(x) dexp(x, rate=1,log=TRUE)
logprior.pi=function(x) dbeta(x,1,1,log=TRUE)
tune=TRUE
max.tune=2e4
file.out="mcmc.txt"
quiet=FALSE



##################################################
## auxiliry functions for gamma parametrisation ##
##################################################
## with shape 'alpha' and rate 'beta'
## E(X) = alpha/beta
## V(X) = alpha/beta^2
## conversely:
## shape: alpha = E(X)^2 / V(X)
## rate: beta = E(X) / V(X)

.gammaToShapeRate <- function(E,V){
    return(c(alpha=E^2/V, beta=E/V))
}

.gammaToMeanVar <- function(alpha,beta){
    return(c(E=alpha/beta, V=alpha/beta^2))
}



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
epidemicMCMC <- function(onset, patch, w, D.patches=NULL, spa.kernel=dexp,
                         n.iter=2e4, sample.every=200,
                         move.R=TRUE, sd.R=0.01, R.ini=1.1,
                         move.delta=TRUE, sd.delta=0.1, delta.ini=1,
                         move.pi=FALSE, sd.pi=0.005, pi.ini=0.5,
                         move.N=TRUE, N.ini=NULL, sd.N=1, move.N.every=50,
                         logprior.delta=function(x) dexp(x, rate=1,log=TRUE),
                         logprior.R=function(x) dgamma(R, shape=2, rate=1, log=TRUE),
                         logprior.pi=function(x) dbeta(x,1,1,log=TRUE),
                         tune=TRUE, max.tune=1e4,
                         file.out="mcmc.txt", quiet=FALSE){

    ## CHECKS / HANDLE ARGUMENTS ##
    x <- data.frame(onset=onset, patch=patch)
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

    ## force matrix
    if(!is.matrix(D.patches)) D.patches <- as.matrix(D.patches)

    ## force null diagonal
    diag(D.patches) <- 0

    ## make sure ordering of patches is consistent
    D.patches <- D.patches[patches,patches]

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
    incid.tot <- sum(incid.vec)


    ## COMPUTE INFECTIOUSNESS (CONSTANT PART OF IT) ##
    ## AUXILIARY FUNCTION:
    ## ## gets infectiousness for day 'date' and
    ## ## a vector of onset dates 'onsets'
    ## get.infec <- function(date, onsets){
    ##     temp <- as.integer(date-onsets)
    ##     out <- sum(w[temp[temp>0]])
    ##     return(out)
    ## }


    ## GET INFECTIOUSNESS ##

    ## GET BETA FOR ALL PATCHES AT TIME T
    ## beta_t_i: infectiousness due to individuals in patch 'i' at time 't'
    ## beta_t^i = sum_{s=1}^{t-1} N_s^i w(t-s)
    ## N is a matrix of incidence (row=dates,col=patches)
    get.betas.t <- function(N,t) return(apply(incid.mat[1:(t-1),,drop=FALSE] * w[(t-1):1], 2, sum))

    ## GET ALL BETAS (excluding t=1)
    get.betas <- function(N) {
        return(t(sapply(2:nrow(N), function(t) get.betas.t(N,t))))
    }



    ## GET LAMBDA FOR ALL PATCHES ##
    ## (force of infection experienced by the patches)
    ## computed using matrix product
    ## (beta %*% P) with:
    ## - beta: infectiousness for day (row) x patches (col)
    ## - P: connectivity from patches (row) to patches (col),
    ## standardized by column
    ## here, beta = R * infect.mat
    ## and P is given by 'spa.kernel(D.pacthes, delta)'
    ## but needs to be column-standardized
    ##
    get.lambdas <- function(N, R, delta){
        return(
            (R* get.betas(N) %*%
             prop.table(spa.kernel(D.patches, delta),2)
             )
            )
    }

    ## faster alternative with pre-computed betas
    get.lambdas.fast <- function(BETAS, R, delta){
        return(
            (R* BETAS %*%
             prop.table(spa.kernel(D.patches, delta),2)
             )
            )
    }


    ## LOG-LIKELIHOOD AND PRIORS ##
    ## all arguments are vectors, except for:
    ## - N (matrix)

    ## GENERAL LIKELIHOOD
    ## p(I|N, pi) p(N|R, delta)
    LL.all <- function(N,pi,R,delta,BETAS){
        return(LL.I(N,pi) + LL.N(N,R,delta,BETAS))
    }

    ## PROBA OF OBSERVED INCIDENCE
    ## p(I|N, pi)
    LL.I <- function(N,pi){
        ## temp <- sum(dbinom(incid.vec, size=as.vector(N), prob=pi, log=TRUE))
        ## if(any(is.na(temp))) {
        ##     cat("\nNA in LL.I")
        ##     print(incid.vec)
        ##     print(as.vector(N))
        ##     print(pi)
        ## }
        return(sum(dbinom(incid.vec, size=as.vector(N), prob=pi, log=TRUE)))
    }

    ## PROBA OF ACTUAL (AUGMENTED) INCIDENCE
    ## p(N | R, delta)
    LL.N <- function(N, R, delta, BETAS){
        rates <- as.vector(get.lambdas.fast(BETAS, R, delta))
        return(sum(dpois(as.vector(N), rates, log=TRUE)))
    }



    ## PROBA OF OSERVED
    ## PRIORS ##
    ## all priors
    logprior.all <- function(delta, R, pi){
        return(logprior.delta(delta) + logprior.R(R) + logprior.pi(pi))
    }



    ## PARAMETER MOVEMENTS ##
    ## MOVE R
    ## movements impact:
    ## p(N | R, delta) p(R)
    R.ACC <- 0
    R.REJ <- 0
    R.move <- function(N, R, delta, BETAS){
        ## generate proposals ##
        newR <- rnorm(n=1, mean=R, sd=sd.R)

        if(newR>=0){
            if(log(runif(1)) <=  (LL.N(N, newR, delta, BETAS) + logprior.R(newR) -
                                  LL.N(N, R, delta, BETAS) - logprior.R(R)
                                  )){
                ## cat("\n(accepted)")
                R <- newR # accept
                R.ACC <<- R.ACC+1
            } else { # reject
                ## cat("\n(rejected)")
                R.REJ <<- R.REJ+1
            }
        } else { # reject
            ## cat("\n(rejected)")
            R.REJ <<- R.REJ+1
        }

        ## return moved vector
        return(R)
    } # end R.move


    ## MOVE SPATIAL PARAM 'DELTA'
    ## movements impact:
    ## p(N | R, delta) p(delta)
    delta.ACC <- 0
    delta.REJ <- 0
    delta.move <- function(N, R, delta, BETAS){
        ## generate proposals ##
        newdelta <- rnorm(n=1, mean=delta, sd=sd.delta)

        if(newdelta>=0){
            if(log(runif(1)) <=  (LL.N(N, R, newdelta, BETAS) + logprior.delta(newdelta) -
                                  LL.N(N, R, delta, BETAS) - logprior.delta(delta))){
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


    ## MOVE PROPORTION OF REPORTED CASES 'pi'
    ## movements impact:
    ## p(I | N, pi) p(pi)
    pi.ACC <- 0
    pi.REJ <- 0
    pi.move <- function(N, pi){
        ## generate proposals ##
        newpi <- rnorm(n=1, mean=pi, sd=sd.pi)

        if(newpi>=0 && newpi<=1){
            if(log(runif(1)) <=  (LL.I(N, newpi) + logprior.pi(newpi) -
                                  LL.I(N, pi) - logprior.pi(pi))){
                pi <- newpi # accept
                pi.ACC <<- pi.ACC+1
            } else { # reject
                pi.REJ <<- pi.REJ+1
            }
        } else { # reject
            pi.REJ <<- pi.REJ+1
        }

        ## return moved vector
        return(pi)
    } # end pi.move



    ## MOVE PROPORTION OF REPORTED CASES 'pi' AND AUGMENTED DATA 'N'
    ## movements impact:
    ## p(I | N, pi) p(pi)
    pi.ACC <- 0
    pi.REJ <- 0
    pi.N.move <- function(N, pi, LAMBDAS){
        ## generate proposals ##
        newpi <- rnorm(n=1, mean=pi, sd=sd.pi)
        newN <- N.move(N, newpi, LAMBDAS)

        if(newpi>=0 && newpi<=1){
            if(log(runif(1)) <=  (LL.I(newN, newpi) + logprior.pi(newpi) -
                                  LL.I(N, pi) - logprior.pi(pi))){
                pi <- newpi # accept
                N <- newN
                pi.ACC <<- pi.ACC+1
            } else { # reject
                pi.REJ <<- pi.REJ+1
            }
        } else { # reject
            pi.REJ <<- pi.REJ+1
        }

        ## return moved vector
        return(list(pi,N))
    } # end pi.N.move


    ## MOVE TRUE, UNOBSERVED INCIDENCE 'N'
    ## uses a Gibbs sampler
    ## (N_t^j-I_t^j) ~ Poisson((1-pi)*lambda_t^j)
    ## all dates but the first one can be moved
    augdat.size <- (n.dates-1)*n.patches
    N.move <- function(N, pi, LAMBDAS){
        ## compute rates for the Poisson distribution
        rates <- (1-pi)*LAMBDAS

        ## draw data
        N[-1,] <- incid.mat[-1,] + rpois(n=augdat.size, lambda=rates)

        ## return moved matrix
        return(N)
    } # end N.move




    ## TUNING ##
    stop.tune.R <- stop.tune.delta <- stop.tune.pi <- FALSE

    ## TUNING FUNCTION FOR R
    R.tune <- function(sd.R){
        ## jumps too large
        if((R.ACC/(R.ACC + R.REJ))<0.2){
            return(sd.R*runif(1,.5,.95))
        }

        ## jumps too small
        if((R.ACC/(R.ACC + R.REJ))>0.5){
            return(sd.R*runif(1,1.05,1.5))
        }

        ## if reaching here, stop tuning R
        stop.tune.R <<- TRUE
    } # end R.tune


    ## TUNING FUNCTION FOR DELTA
    delta.tune <- function(sd.delta){
        ## jumps too large
        if((delta.ACC/(delta.ACC + delta.REJ))<0.2){
            return(sd.delta*runif(1,.5,.95))
        }

        ## jumps too small
        if((delta.ACC/(delta.ACC + delta.REJ))>0.5){
            return(sd.delta*runif(1,1.05,1.5))
        }

        ## if reaching here, stop tuning delta
        stop.tune.delta <<- TRUE
    } # end delta.tune


    ## TUNING FUNCTION FOR PI
    pi.tune <- function(sd.pi){
        ## jumps too large
        if((pi.ACC/(pi.ACC + pi.REJ))<0.2){
            return(sd.pi*runif(1,.5,.95))
        }

        ## jumps too small
        if((pi.ACC/(pi.ACC + pi.REJ))>0.5){
            return(sd.pi*runif(1,1.05,1.5))
        }

        ## if reaching here, stop tuning pi
        stop.tune.pi <<- TRUE
    } # end pi.tune


    ## INITIALIZE VALUES
    R <- R.ini
    delta <- delta.ini
    pi <- pi.ini
    N <- round(incid.mat/pi)
    BETAS <- get.betas(N)
    LAMBDAS <- get.lambdas.fast(BETAS, R, delta)

    ## check if initial log likelihood is -Inf
    if(!is.finite(LL.all(N,pi,R,delta,BETAS))){
        cat("\nProblem: initial log-likelihood is -Inf")
        cat("\nDetails of the log-likelihood components")
        cat("\np(I|N, pi): ", LL.I(N, pi))
        cat("\np(N|R, delta): ", LL.N(N, R, delta, BETAS))
        cat("\n")
    }

    ## check if initial log prior is -Inf
    if(!is.finite(logprior.all(delta, R, pi))){
        cat("\nProblem: initial log-prior is -Inf")
        cat("\nDetails of the log-prior components")
        cat("\np(delta): ", logprior.pi(delta))
        cat("\np(R): ", logprior.R(R))
        cat("\np(pi): ", logprior.pi(pi))
        cat("\n")
    }

    ## basic message
    if(!quiet && tune) cat("Starting tuning proposal distributions...\n")

    COUNTER <- 0
    KEEPTUNING <- tune

    while(KEEPTUNING && (COUNTER<=max.tune)){
        ## update counter
        COUNTER <- COUNTER + 1

        ## move stuff ##
        ## move R if needed
        if(move.R) R <- R.move(N, R, delta, BETAS)

        ## move delta if needed
        if(move.delta) delta <- delta.move(N, R, delta, BETAS)

        ## move pi if needed
        if(move.pi) {
            temp <- pi.N.move(N, pi, LAMBDAS)
            pi <- temp[[1]]
            N <- temp[[2]]
        }

        ## move N if needed
        if(move.N && (COUNTER %% move.N.every < 1)) {
            N <- N.move(N, pi, LAMBDAS)
            BETAS <- get.betas(N)
        }

        ## change parameters of proposal distributions ##
        ## (every 100 iterations)
        if(COUNTER %% 100 ==0){
            ## update SDs
            if(move.R) sd.R <- R.tune(sd.R) else stop.tune.R <- TRUE
            if(move.delta) sd.delta <- delta.tune(sd.delta) else stop.tune.delta <- TRUE
            if(move.pi) sd.pi <- pi.tune(sd.pi) else stop.tune.pi <- TRUE

            ## check if we should stop
            KEEPTUNING <- !(all(stop.tune.R, stop.tune.delta, stop.tune.pi))
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

        ## acceptance rates for pi
        if(move.pi){
            cat("\ntuned acceptance rate for pi: ", pi.ACC/(pi.ACC+pi.REJ))
            cat("\naccepted: ", pi.ACC)
            cat("\nreject: ", pi.REJ)
            cat("\ntotal: ", pi.ACC + pi.REJ)
            cat("\n")
        }

    }



    ## MCMC ##
    ## BASIC HEADER
    header <- "step\tpost\tlikelihood\tprior\tR\tdelta\tpi"
    cat(header, file=file.out)

    ## add first line
    ## temp: c(loglike, logprior)
    temp <- c(LL.all(N,pi,R,delta,BETAS), logprior.all(delta=delta, R=R, pi=pi))

    ## check that initial LL is not -Inf
    if(!is.finite(temp[1])) warning("Initial likelihood is zero")

    ## write to file
    cat("\n", file=file.out, append=TRUE)
    cat(c(1, sum(temp), temp, R, delta, pi), sep="\t", append=TRUE, file=file.out)

    ## basic message
    if(!quiet) cat("\nStarting MCMC: 1")

    ## mcmc ##
    for(i in 1:n.iter){
        ## move stuff ##
        ## move R if needed
        if(move.R) R <- R.move(N, R, delta, BETAS)

        ## move delta if needed
        if(move.delta) delta <- delta.move(N, R, delta, BETAS)

        ## move pi if needed
        if(move.pi) {
            temp <- pi.N.move(N, pi, LAMBDAS)
            pi <- temp[[1]]
            N <- temp[[2]]
        }

        ## move N if needed
        if(move.N && (i %% move.N.every < 1)) {
            N <- N.move(N, pi, LAMBDAS)
            BETAS <- get.betas(N)
        }

        ## if retain this sample ##
        if(i %% sample.every ==0){
            temp <- c(LL.all(N,pi,R,delta,BETAS), logprior.all(delta=delta, R=R, pi=pi))
            cat("\n", file=file.out, append=TRUE)

            ## add posterior
            cat(c(i, sum(temp), temp, R, delta, pi), sep="\t", append=TRUE, file=file.out)

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
    out$AR <- data.frame(accept=c(R.ACC, delta.ACC, pi.ACC),
                         reject=c(R.REJ, delta.REJ, pi.REJ))
    out$AR$prop <- out$AR$accept / (out$AR$accept + out$AR$reject)
    out$AR$sd <- c(sd.R, sd.delta, sd.pi)
    rownames(out$AR) <- c("R","delta","pi")

    ## add class, and return ##
    class(out) <- c("list", "epidemicMCMC")
    return(out)

} # end epidemicMCMC





