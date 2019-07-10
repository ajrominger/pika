#' @title Log likelihood methods for objects of class \code{sad}
#'  
#' @description Compute log likelihood and likelihood-based z-scores for objects of class \code{sad}.
#' 
#' @details \code{logLikZ.sad} 
#' simulates from a fitted SAD and calculates the likelihood of these simulated data sets. 
#' The distribution of these values is compared against the likelihood of the data to obtain 
#' a z-score.
#' 
#' @param x an object of class \code{sad}
#' @param ... arguments to be passed to methods
#' 
#' @keywords species abundance, SAD, log likelihood
#' @rdname logLik
#' @export
#' 
#' @examples
#' x <- rtpois(1000, 2)
#' sad1 <- fitSAD(x, models='tpois')
#' logLik(sad1[[1]])
#' 
#' @return \code{logLik} returns an object of class \code{logLik} with attributes "\code{df}" and "\code{nobs}". 
#' \code{logLikZ} returns a list with elements
#' \describe{
#'    \item{z}{The z-score}
#'    \item{obs}{The observed log-likelihood of the actual data set}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD, logLik
# @references

logLik.sad <- function(x) {
    out <- x$ll
    attr(out, 'df') <- x$df
    attr(out, 'nobs') <- x$nobs
    class(out) <- 'logLik'
    
    return(out)
}


#======================================================
#' @rdname logLik
#' @export 

logLikZ <- function(x, ...) {
    UseMethod('logLikZ')
}

#' @rdname logLik
#' @export 
#' @importFrom stats logLik sd

logLikZ.sad <- function(x, ...) {
    lik.obs <- logLik(x)
    n <- x$nobs
    
    ## hypothetical distribution of probabilities
    p0 <- .p0(x)
    
    ## hypothetical mean and var
    m <- sum(p0 * exp(p0)) * n
    v <- sum((m/n - p0)^2 * exp(p0)) * n
    
    ## z^2-value
    z <- ((lik.obs - m) / sqrt(v))^2
    
    return(as.numeric(z))
}


## helper function to compute distribution of potential probabilities under the 
## hypothesis that the SAD of interest is true

.p0 <- function(x) {
    n0 <- 1:10^5
    dfun <- getdfun(x)
    
    p0 <- dfun(n0, log = TRUE)
    p0 <- p0[is.finite(p0)]
    
    if(exp(p0[length(p0)]) > .Machine$double.eps^0.75) {
        n0add <- (n0[length(p0)] + 1):10^6
        p0add <- dfun(n0add, log = TRUE)
        p0add <- p0add[is.finite(p0add)]
        
        p0 <- c(p0, p0add)
    }
    
    return(p0)
}
