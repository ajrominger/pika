#' @title Fit a model for a species abundance distribution
#'  
#' @description \code{fitSAD} uses maximum likelihood to fit one or several SAD models to data
#' 
#' @details The default behavior fits all models availible in \code{pika}:
#' \describe{
#'  \item{\code{lseries}}{Log series}
#'  \item{\code{plnorm}}{Poisson log normal}
#'  \item{\code{stick}}{Broken stick}
#'  \item{\code{tnegb}}{Truncated negative binomial}
#'  \item{\code{tpois}}{Truncated Poisson}
#' }
#' Exact solutions are availible for \code{tpois} and \code{stick}. In the case of \code{lseries}, 
#' optimization is used to solve \eqn{\bar{x} = \frac{-1}{log(1-e^{-\beta})}\frac{e^{-\beta}}{1-e^{-\beta}}}. 
#' For \code{plnorm} and \code{tnegb} optimization is used to maximize the log likelihood functions based 
#' directly on \code{dplnorm} and \code{dtnegb}, respectively.
#' 
#' @param x vector of integers representing a sample of species abundances
#' @param models character vector of models to be fit (can be more than one model, see examples)
#' @param keepData logical, should the data be stored with the fitted \code{sad} object
#' 
#' @keywords Maximum likelihood, species abundance, SAD
#' @export
#' 
#' @examples
#' x <- rtpois(1000, 2)
#' fitSAD(x, models='tpois')
#' 
#' @return A list containing named objects of class \code{sad} with elements
#' \describe{
#'   \item{\code{MLE}}{The maximum likelihood estimate(s)}
#'   \item{\code{ll}}{The maximized log likelihood}
#'   \item{\code{df}}{The associated degrees of freedom}
#'   \item{\code{nobs}}{The associated number of observations}
#'   \item{\code{model}}{The name of the model fit}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso logLik.sad, dlseries, dplnorm, dstick, dtnegb, dtpois
# @references 

fitSAD <- function(x, models=c('lseries', 'plnorm', 'stick', 'tnegb', 'tpois'), keepData=FALSE) {
    if(length(models) > 1) {
        out <- lapply(models, function(m) fitSAD(x, m, keepData)[[1]])
        # return(out)
    } else {
        out <- switch(models,
                      'lseries' = .fitlseries(x),
                      'plnorm' = .fitPlnorm(x),
                      'stick' = .fitStick(x),
                      'tnegb' = .fitTnegb(x),
                      'tpois' = .fitTpois(x))
        
        out$model <- models
        if(keepData) out$data <- x
        
        class(out) <- 'sad'
    }
    
    if(length(models) == 1) {
        out <- list(out)
    }
    
    names(out) <- models
    
    return(out)
}



## ========================================================
## helper functions to return likelihood estimates for SADs
## ========================================================

## MLE for log series
#' @export
.fitlseries <- function(x) {
    xbar <- mean(x)
    fun <- function(b) -1/log(1 - exp(-b)) * exp(-b)/(1 - exp(-b)) - xbar
    fit <- uniroot(fun, .lseriesSolInt(xbar), tol=.Machine$double.eps)
    return(list(MLE=fit$root, ll=sum(dlseries(x, fit$root, log=TRUE)), df=1, nobs=length(x)))
}

## function to find bounding interval of MLE for logseries
#' @export
.lseriesSolInt <- function(x) {
    up <- 1.1*exp(2*(x-0.1)^-6) * (x-0.1)^-1.25
    lo <- 0.5*exp(2*(x+0.1)^-5) * (x+0.1)^-1.5
    return(c(lo, up))
}

## MLE for Poisson log normal
#' @export
.fitPlnorm <- function(x) {
    fun <- function(par) {
        out <- -sum(dplnorm(x, par[1], par[2], log=TRUE))
        if(!is.finite(out)) out <- .Machine$double.xmax^0.25
        return(out)
    }
    
    init.mu <- mean(log(x))
    init.sig <- sd(log(x))
    
    fit <- optim(c(init.mu, init.sig), fun, method='L-BFGS-B', 
                 lower=c(-10*init.mu, .Machine$double.eps), upper=2*c(init.mu, init.sig))
    
    return(list(MLE=fit$par, ll=-fit$value, df=2, nobs=length(x)))
}


## MLE for broken stick
#' @export
.fitStick <- function(x) {
    mle <- 1/mean(x)
    
    return(list(MLE=mle, ll=sum(dstick(x, mle, log=TRUE)), df=1, nobs=length(x)))
}

## MLE for truncated negative binomial
#' @export
.fitTnegb <- function(x) {
    fun <- function(par) {
        par[2] <- exp(par[2])
        out <- -sum(dtnegb(x, par[1], par[2], log=TRUE))
        if(!is.finite(out)) out <- .Machine$double.xmax^0.25
        
        return(out)
    }
    
    init.mu <- mean(x)
    if(var(x) < mean(x)) {
        init.k <- log(10000)
    } else {
        init.k <- log(mean(x)^2 / (var(x) - mean(x)))
    }
    
    fit <- optim(c(init.mu/2, init.k/2), fun, method='L-BFGS-B', 
                 lower=c(.Machine$double.eps, log(.Machine$double.eps)), upper=c(init.mu*10, init.k+10))
    fit$par[2] <- exp(fit$par[2])
    
    return(list(MLE=fit$par, ll=-fit$value, df=2, nobs=length(x)))
}


## MLE for truncated Poisson
#' @export
.fitTpois <- function(x) {
    mle <- mean(x) + gsl::lambert_W0(-mean(x)*exp(-mean(x)))
    
    return(list(MLE=mle, ll=sum(dtpois(x, mle, log=TRUE)), df=1, nobs=length(x)))
}

