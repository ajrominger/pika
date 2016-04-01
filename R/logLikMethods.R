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
#' @seealso fitSAD
# @references

logLik.sad <- function(x) {
    out <- x$ll
    attr(out, 'df') <- x$df
    attr(out, 'nobs') <- x$nobs
    class(out) <- 'logLik'
    
    return(out)
}


#==============================================================================
#' @rdname logLik
#' @export 

logLikZ <- function(x, ...) {
    UseMethod('logLikZ')
}

#' @rdname logLik
#' @export 

#' @importFrom stats logLik sd
logLikZ.meteDist <- function(x, nrep=999, return.sim=FALSE, ...) {
    lik.obs <- logLik(x)
    state.var <- sum(x$data)
    
    lik.sim <- c()
    cat('simulating data that conform to state variables: \n')
    for(i in 1:10) {
        cat(sprintf('attempt %s \n', i))
        this.sim <- replicate(100*nrep, {
            new.dat <- x$r(length(x$data))
            if(abs(sum(new.dat) - state.var) < 0.001*state.var) {
                return(NA)
            } else {
                return(sum(x$d(new.dat, log=TRUE)))
            }
        })
        
        lik.sim <- c(lik.sim, this.sim[!is.na(this.sim)])
        if(length(lik.sim) >= nrep) break
    }
    
    if(length(lik.sim) >= nrep) {
        lik.sim <- c(lik.sim[1:nrep], lik.obs)
    } else {
        warning(sprintf('%s (not %s as desired) simulated replicates found that match the state variables', 
                        length(lik.sim), nrep))
        lik.sim <- c(lik.sim, lik.obs)
    }
    
    z <- (lik.obs-mean(lik.sim))/sd(lik.sim)
    
    if(!return.sim) lik.sim <- NULL
    
    return(list(z=z, 
                obs=lik.obs,
                sim=lik.sim))
    
}
