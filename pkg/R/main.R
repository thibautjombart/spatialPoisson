

## to simulate data
x <- data.frame(onset=sample(as.Date("2014-01-01")+0:10, 30, replace=TRUE), patch=sample(c('a','b','c'), replace=TRUE, 30))




model <- function(x, D.patches=NULL, w){
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
    ## force diagonal to zero
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

    ## vectorized format for faster computations in dpois
    incid.vec <- as.vector(incid[-1,]) # first date is removed (infectiousness is 0)


    ## FUNCTIONS TO GET RATES AND LOG-LIKELIHOOD ##
    ## POISSON LIKELIHOOD
    ## p(y|rates)
    get.LL <- function(y, rates){
        return(sum(dpois(y, rates, log=TRUE)))
    }

    ## GET BETA FOR A PATCH ##
    ## get infectiousness for day 'date' and
    ## a vector of onset dates 'onsets'
    get.infec <- function(date, onsets){
        temp <- as.integer(date-onsets)
        out <- sum(w[temp[temp>0]])
        return(out)
    }

    ## 'infec.mat' contains the sum of infectiousness for each day (row) and patch (column)
    ## basic matrix
    infec.mat <- t(sapply(all.dates, function(t) tapply(x$onset, x$patch, function(e) get.infec(t, e))))

    ## vectorized format for faster computations in dpois
    infec.vec <- as.vector(infec.mat[-1,]) # first date is removed (infectiousness is 0)


}
