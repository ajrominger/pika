#' @title Compute residuals and mean squared error between observed and predicted SADs.
#'
#' @description
#' \code{residuals.sad} computes residuals between predicted and observed SAD, either the rank abundance or CDF. 
#' \code{mse.sad} computes mean squared error for rank or cdf between predicted and observed SAD.
#' \code{mseZ.sad} computes the z-score of mean squared error based on parametric bootstrapping.
#'
#' @details
#' Residuals are calculated either for the rank abundances or the cumulative distribution function. 
#' This can be done on a log scale and/or as the proportional error.  The mean squared error is then 
#' simply the mean of squared residuals of either the rank or cumulative distribution function. 
#' z-scores are caculated by simulating from a fitted SAD and calculating the MSE between the simulated 
#' data sets and the SAD prediction. The distribution of these values is compared against the MSE of the 
#' data to obtain a z-score. 
#' 
#' @param x an \code{sad} object
#' @param type 'rank' or 'cumulative'
#' @param relative logical, if true use relative MSE; defaults to TRUE
#' @param log logical, if TRUE calculate MSE on log scale; defaults to FALSE
#' @param nrep number of simulations from the fitted METE distribution 
#' @param return.sim logical; return the simulated liklihood values
#' @param ... arguments to be passed to methods
#' 
#' @rdname residuals
#' @export
#' 
#' @examples
#' 'stub'
#' 
#' @return 
#' \describe{
#'  \item{\code{residuals}}{returns a numeic vector giving residuals for each data point; 
#'  for \code{type='rank'} this will be a vector of length equal to the number of species 
#'  but for \code{type='cumulative'} this will be a vector of length equal to the number 
#'  of unique abundances.}
#'  \item{\code{mse}}{returns a numeric vector of length one giving the value of the mean 
#'  squared error.}
#'  \item{\code{mseZ}}{returns a list with elements: \code{z}, the z-score; \code{obs}, 
#'  the observed mean squared error; \code{sim}, \code{nrep} simulated values (\code{NULL} 
#'  if \code{return.sim=FALSE}.)}
#' }
#' 
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso

residuals.sad <- function(x, type=c("rank","cumulative"),
                          relative=TRUE, log=FALSE, ...) {
    type <- match.arg(type, choices=c("rank","cumulative"))
    
    if(type=="rank") {
        obs <- sort(x$data, decreasing=TRUE)
        pred <- sad2Rank(x)
        if(log) {
            obs <- log(obs)
            pred <- log(pred)
        }
    } else if(type=="cumulative") {
        obs <- .ecdf(x$data)
        pfun <- getpfun(x)
        
        pred <- pfun(obs[, 1], log=log)
        if(log) obs <- log(obs[, 2])
        else obs <- obs[,2]
    }
    
    out <- obs - pred
    if(relative) out <- out/abs(pred)
    
    return(out)
}


## ================================================

#' @rdname residuals
#' @export

mse <- function(x, ...) {
    UseMethod('mse')
}

#' @rdname residuals
#' @export 
#' @importFrom stats residuals
mse.sad <- function(x, ...) {
    resid <- residuals(x, ...)
    
    return(mean(resid^2))
}


## ================================================

#' @rdname residuals
#' @export

mseZ <- function(x, ...) {
  UseMethod('mseZ')
}

#' @rdname residuals
#' @export 
#' @importFrom stats sd
mseZ.sad <- function(x, nrep, return.sim=FALSE, ...) {
    this.mseFun <- function(x) mse(x, ...)
    mse.obs <- this.mseFun(x)
    
    rfun <- getrfun(x)
    n <- x$nobs
    newx <- x
    
    mse.sim <- replicate(nrep, {
        newx$data <- rfun(n)
        this.mseFun(newx)
    })
    
    z <- ((mse.obs - mean(mse.sim)) / sd(mse.sim))^2
    
    if(return.sim) {
        mse.sim <- ((mse.sim - mean(mse.sim)) / sd(mse.sim))^2
    } else {
        mse.sim <- NULL
    }
    
    return(list(z=z, obs=mse.obs, sim=mse.sim))
}
