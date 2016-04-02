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
#' @param nrep number of simulations from the fitted METE distribution 
#' @param return.sim logical; return the simulated liklihood values
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
#'    \item{sim}{\code{nrep} simulated values if \code{return.sim=TRUE}, \code{NULL} otherwise.}
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

logLikZ.sad <- function(x, nrep=1000, return.sim=FALSE, ...) {
    lik.obs <- logLik(x)
    n <- x$nobs
    rfun <- getrfun(x)
    dfun <- getdfun(x)
    
    lik.sim <- replicate(nrep, {
        newx <- rfun(n)
        sum(dfun(newx, log=TRUE))
    })
    
    z <- ((lik.obs - mean(lik.sim)) / sd(lik.sim))^2
    
    if(return.sim) {
        lik.sim <- ((lik.sim - mean(lik.sim)) / sd(lik.sim))^2
    } else {
        lik.sim <- NULL
    }
    
    return(list(z=z, obs=lik.obs, sim=lik.sim))
}
