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
#' @param relative logical; if true use relative MSE
#' @param log logical; if TRUE calculate MSE on logged distirbution. If FALSE use arithmetic scale.
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
        obs <- x$data
        pred <- sad2Rank(x)
    } else if(type=="cumulative") {
        obs <- .ecdf(x$data)
        
        pfun <- switch(x$model, 
                       'fish' = function(q) pfish(q, x$MLE, log=log),
                       'plnorm' = function(q) pplnorm(q, x$MLE[1], x$MLE[2], log=log),
                       'stick' = function(q) pstick(q, x$MLE, log=log),
                       'tnegb' = function(q) ptnegb(q, x$MLE[1], x$MLE[2], log=log),
                       'tpois' = function(q) ptpois(q, x$MLE, log=log))
        
        pred <- pfun(obs[, 1])
        if(log) obs <- log(obs[, 2])
        else obs <- obs[,2]
    }
    
    out <- obs - pred
    if(relative) out <- out/abs(pred)
    
    return(out)
}


#=======================================================================

#' @rdname residuals
#' @export

mse <- function(x, ...) {
    UseMethod('mse')
}

# @rdname residuals
# @export 
# @importFrom stats residuals 
# mse.meteDist <- function(x, type=c("rank","cumulative"),
#                          relative=TRUE, log=FALSE, ...) {
#     type <- match.arg(type, choices=c("rank","cumulative"))
#     
#     resid <- residuals(x, type, relative, log)
#     
#     return(mean(resid^2))
# }

#' @rdname residuals
#' @export 
#' @importFrom stats residuals
mse.sad <- function(x,...) {
    resid <- residuals(x, ...)
    
    return(mean(resid^2))
}


#=======================================================================
#' @rdname residuals
#' @export

mseZ <- function(x, ...) {
  UseMethod('mseZ')
}

#' @rdname residuals
#' @export 
#' @importFrom stats sd
mseZ.sad <- function(x, nrep, return.sim=FALSE,
                     type=c("rank","cumulative"), 
                     relative=TRUE, log=FALSE, ...) {
    
    if(type=='rank') {
        rad <- meteDist2Rank(x)
        thr <- function(dat) {
            res <- sort(dat, TRUE) - rad
            if(relative) res <- res/rad
            
            return(mean(res^2))
        }
    } else {
        obs <- .ecdf(dat)
        thr <- function(dat) {
            if(log) obs[, 2] <- log(obs[, 2])
            pred <- x$p(obs[, 1], log.p=log)
            res <- obs[, 2] - pred
            if(relative) res <- res/abs(pred)
            
            return(mean(res^2))
        }
    }
    
    mse.obs <- mse.meteDist(x, type, relative, log)
    state.var <- sum(x$data)
    
    mse.sim <- c()
    cat('simulating data that conform to state variables: \n')
    for(i in 1:10) {
        cat(sprintf('attempt %s \n', i))
        this.sim <- replicate(100*nrep, {
            new.dat <- x$r(length(x$data))
            if(abs(sum(new.dat) - state.var) < 0.001*state.var) {
                return(NA)
            } else {
                return(thr(new.dat))
            }
        })
        
        mse.sim <- c(mse.sim, this.sim[!is.na(this.sim)])
        if(length(mse.sim) >= nrep) break
    }
    
    if(length(mse.sim) >= nrep) {
        mse.sim <- c(mse.sim[1:nrep], mse.obs)
    } else {
        warning(sprintf('%s (not %s as desired) simulated replicates found that match the state variables', 
                        length(lik.sim), nrep))
        lik.sim <- c(mse.sim, mse.obs)
    }
    
    if(return.sim) {
        return(list(z=(mse.obs-mean(mse.sim))/sd(mse.sim), 
                    obs=mse.obs,
                    sim=mse.sim))
    } else {
        return(list(z=(mse.obs-mean(mse.sim))/sd(mse.sim)))
    }
}
