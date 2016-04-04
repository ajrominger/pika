#' @title Broken stick
#'  
#' @description \code{dstick} gives the probability mass function, \code{pstick} gives the cumulative mass 
#' function, \code{qstick} the quantile function, \code{rstick} randome number generation
#' 
#' @details These functions assume infinite support of the Broken stick model from [1, Inf). Parameter \code{r} corresponds roughly to the probability of dividing the resource, but really should be interpreted in terms of the underlying PMF: \eqn{(1-r)^{x-1}r}
#' 
#' @param x vector of integers for which to return the probability
#' @param q vector of integers for which to return the cumulative probability
#' @param p vector of probabilities for which to return the quantile
#' @param n number of random replicates
#' @param r parameter of the broken stick model, see Details
#' @param log logical, should the log probability be used
#' @param lower.tail logical, should the lower tail be used
#' 
#' @keywords Broken stick, geometric distribution, species abundance, SAD
#' @export
#' 
#' @examples
#' 
#' dstick(1:10, 0.5)
#' 
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD, dgeom
# @references 

#' @rdname Stick

dstick <- function(x, r, log=FALSE) {
    out <- (1-r)^(x-1) * r
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer x = %s', bad))
        }
        
        out[x %% 1 != 0] <- 0
    }
    
    out[x < 1] <- 0
    
    if(log) out <- log(out)
    return(out)
}


#' @export
#' @rdname Stick

pstick <- function(q, r, lower.tail=TRUE, log=FALSE) {
    out <- 1 - (1-r)^q
    
    if(any(q %% 1 != 0)) {
        for(bad in q[q %% 1 != 0]) {
            warning(sprintf('non-integer q = %s', bad))
        }
        
        out[q %% 1 != 0] <- 0
    }
    
    out[q < 1] <- 0
    
    if(!lower.tail) out <- 1 - out
    if(log) out <- log(out)
    return(out)
}


#' @export
#' @rdname Stick

qstick <- function(p, r, lower.tail=TRUE, log=FALSE) {
    if(log) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    
    out <- .stickcdfinv(p, r)
    
    if(any(is.nan(out))) {
        warning('NaNs produced')
    }
    
    return(out)
}


#' @export
#' @rdname Stick

rstick <- function(n, r) {
    return(rgeom(n, r) + 1)
}


## =================================
## helper functions
## =================================

#' @export
## inverse cdf of the broken stick
.stickcdfinv <- function(p, r) {
    log(1 - p) / log(1 - r)
}
