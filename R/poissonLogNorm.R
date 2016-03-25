#' @title Poisson log-normal
#'  
#' @description \code{dplnorm} gives the probability mass function, \code{pplnorm} gives the cumulative mass 
#' function, \code{qplnorm} the quantile function, \code{rplnorm} randome number generation
#' 
#' @details These functions assume infinite support of the Poisson log-normal from [1, Inf).
#' 
#' @param x vector of integers for which to return the probability
#' @param q vector of integers for which to return the cumulative probability
#' @param p vector of probabilities for which to return the quantile
#' @param n number of random replicates
#' @param mu mean log abundance, see Details
#' @param sig standard deviation of log abundance, see Details
#' @param log logical, should the log probability be used
#' @param lower.tail logical, should the lower tail be used
#' 
#' @keywords Poisson log normal, species abundance, SAD
#' @export
#' 
#' @examples
#' 
#' dplnorm(1:10, 0.5, 0.1)
#' 
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso 
# @references 

#' @rdname PoisLogNormal

dplnorm <- function(x, mu, sig, log=FALSE) {
    out <- poilog::dpoilog(x, mu, sig) / (1 - poilog::dpoilog(0, mu, sig))
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer x = %s', bad))
        }
        
        out[x %% 1 != 0] <- 0
    }
    
    if(log) out <- log(out)
    return(out)
}


#' @rdname PoisLogNormal

pplnorm <- function(q, mu, sig, lower.tail=TRUE, log=FALSE) {
    out <- 1 - (1-r)^q
    
    if(any(q %% 1 != 0)) {
        for(bad in q[q %% 1 != 0]) {
            warning(sprintf('non-integer q = %s', bad))
        }
        
        out[q %% 1 != 0] <- 0
    }
    
    if(!lower.tail) out <- 1 - out
    if(log) out <- log(out)
    return(out)
}

#' @rdname PoisLogNormal

qplnorm <- function(p, r, lower.tail=TRUE, log=FALSE) {
    if(log) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    
    out <- .stickcdfinv(p, r)
    
    if(any(is.nan(out))) {
        warning('NaNs produced')
    }
    
    return(out)
}


#' @rdname PoisLogNormal

rplnorm <- function(n, r) {
    lat <- rlnorm()
}


## =================================
## helper functions
## =================================

## inverse cdf of the broken stick
.stickcdfinv <- function(p, r) {
    log(1 - p) / log(1 - r)
}
