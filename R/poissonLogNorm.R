#' @title Poisson log-normal
#'  
#' @description \code{dplnorm} gives the probability mass function, \code{pplnorm} gives the cumulative mass 
#' function, \code{qplnorm} the quantile function, \code{rplnorm} randome number generation
#' 
#' @details These functions assume infinite support of the Poisson log-normal from [1, Inf).  Calculations are based on code from the \code{poilog} package, see References
#' 
#' @param x vector of integers for which to return the probability
#' @param q vector of integers for which to return the cumulative probability
#' @param p vector of probabilities for which to return the quantile
#' @param n number of random replicates
#' @param mu mean log abundance
#' @param sig standard deviation of log abundance
#' @param log logical, should the log probability be used
#' @param lower.tail logical, should the lower tail be used
#' 
#' @keywords Poisson log normal, species abundance, SAD
#' @export
#' 
#' @examples
#' dplnorm(1:10, 0.5, 0.1)
#' 
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso dpois, dlnorm, dtpois
#' @references Engen, S., R. Lande, T. Walla & P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional Poisson lognormal species abundance model. American Naturalist 160: 60-73.

#' @rdname PoisLogNormal

dplnorm <- function(x, mu, sig, log=FALSE) {
    out <- numeric(length(x))
    out[x %% 1 == 0 & x >= 1] <- poilog::dpoilog(x[x %% 1 == 0 & x >= 1], mu, sig) / (1 - poilog::dpoilog(0, mu, sig))
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer x = %s', bad))
        }
    }
    
    if(log) out <- log(out)
    return(out)
}


#' @export
#' @rdname PoisLogNormal

pplnorm <- function(q, mu, sig, lower.tail=TRUE, log=FALSE) {
    out <- numeric(length(q))
    newq <- floor(q)
    newq[newq < 1] <- 1
    
    if(length(q) > 1) {
        temp <- cumsum(poilog::dpoilog(1:max(newq), mu, sig)) / (1 - poilog::dpoilog(0, mu, sig))
        out <- temp[newq]
    } else {
        out <- sum(poilog::dpoilog(1:newq, mu, sig)) / (1 - poilog::dpoilog(0, mu, sig))
    }
    
    out[q < 1] <- 0
    
#     if(any(q %% 1 != 0)) {
#         for(bad in q[q %% 1 != 0]) {
#             warning(sprintf('non-integer q = %s', bad))
#         }
#     }
    
    if(!lower.tail) out <- 1 - out
    if(log) out <- log(out)
    return(out)
}

#' @export
#' @rdname PoisLogNormal

qplnorm <- function(p, mu, sig, lower.tail=TRUE, log=FALSE) {
    if(log) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    
    out <- .plnormcdfinv(p, mu, sig)
    
    if(any(is.nan(out))) {
        warning('NaNs produced')
    }
    
    return(out)
}

#' @export
#' @rdname PoisLogNormal

rplnorm <- function(n, mu, sig) {
    N <- 100 * n / (1 - dpois(0, exp(mu + sig^2/2)))
    lat <- rlnorm(N, mu, sig)
    rel <- rpois(N, lat)
    rel <- rel[rel > 0]
    
    if(length(rel) < n) warning(sprintf('could not find %s unique random variates, using bootstrapping', n))
    
    return(sample(rel, n, rep=ifelse(n < length(rel), FALSE, TRUE)))
}


## =================================
## helper functions
## =================================

## inverse cdf of the poisson log normal
#' @export
.plnormcdfinv <- function(p, mu, sig) {
    this.cdf <- c(0, cumsum(poilog::dpoilog(1:10000, mu, sig)) / (1 - poilog::dpoilog(0, mu, sig)))
    approx(x=this.cdf, y=0:10000,
           xout=p, method='constant', yleft=NaN, yright=Inf, f=1)$y
}
