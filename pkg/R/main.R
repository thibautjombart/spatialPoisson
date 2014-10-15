

## to simulate data
x <- data.frame(onset=sample(as.Date("2014-01-01")+0:10, 30, replace=TRUE), patch=sample(c('a','b','c'), replace=TRUE, 30))




model <- function(x){
    ## CHECKS ##
    ## if(!require(OutbreakTools)) stop("OutbreakTools is needed")
    x <- na.omit(x)
    x$patch <- factor(x$patch)
    patches <- levels(x$patch)
    n.patches <- length(patches)

    ## COMPUTE INCIDENCE ##
    date.range <- range(x$onset)
    all.dates <- seq(date.range[1],date.range[2],by=1)
    n.dates <- length(all.dates)

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
    incid <- data.frame(date=rep(all.dates,n.patches),
                        patch=rep(patches, n.dates),
                        incidence=unlist(temp))


    ## FUNCTION TO GET THE LOG-LIKELIHOOD ##
    ## Poisson likelihood
    ## p(y|rates)
    get.LL <- function(y, rates){
        return(sum(dpois(y, rates, log=TRUE)))
    }

}
