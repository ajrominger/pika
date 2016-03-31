#' @title Fit a model for a species abundance distribution
#'  
#' @description \code{fitSAD} uses maximum likelihood to fit one or several SAD models to data
#' 
#' @details The default behavior fits all models availible in \code{pika}
#' 
#' @param x vector of integers representing a sample of species abundances
#' @param models character vector of models to be fit (can be more than one model, see examples)
#' 
#' @keywords Maximum likelihood, species abundance, SAD
#' @export
#' 
#' @examples
#' 
#' dplnorm(1:10, 0.5, 0.1)
#' 
#' @return An object of class \code{logLik} with attributes "\code{df}" and "\code{nobs}."
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso logLik
# @references 

fitSAD <- function(x, models=c('fish', 'plnorm', 'stick', 'tnegb', 'tpois')) {
    if(length(models) > 1) {
        ## sapply over models, fitting each one...could be done recursively
    } else {
        out <- switch(models,
                      'fish' = .fitFish(x),
                      'plnorm' = .fitPlnorm(x),
                      'stick' = .fitStick(x),
                      'tnegb' = .fitTnegb(x),
                      'tpois' = .fitTpois(x))
        ## extract estimate and df from out, nobs from length(x)
    }
}



## ========================================================
## helper functions to return likelihood estimates for SADs
## ========================================================


## MLE for fisher log series
.fitFish <- function(x) {
    xbar <- mean(x)
    fun <- function(b) -1/log(1 - exp(-b)) * exp(-b)/(1 - exp(-b)) - xbar
    fit <- uniroot(fun, .fishSolInt(xbar), tol=.Machine$double.eps)
    return(list(MLE=fit$root, ll=sum(dfish(x, fit$root, log=TRUE)), df=1, nobs=length(x)))
}

## function to find bounding interval of MLE for fisher
.fishSolInt <- function(x) {
    up <- 1.1*exp(2*(x-0.1)^-6) * (x-0.1)^-1.25
    lo <- 0.5*exp(2*(x+0.1)^-5) * (x+0.1)^-1.5
    return(c(lo, up))
}

## MLE for Poisson log normal 
.fitPlnorm <- function(x) {
    fun <- function(par) {
        par[2] <- exp(par[2])
        return(-sum(dplnorm(x, par[1], par[2], log=TRUE)))
    }
    
    init.mu <- mean(log(x))
    init.sig <- log(sd(log(x)))
    
    fit <- optim(c(init.mu, init.sig), fun, control=list(reltol=.Machine$double.eps))
    fit$par[2] <- exp(fit$par[2])
    
    return(list(MLE=fit$par, ll=-fit$value, df=2, nobs=length(x)))
}


.fitStick <- function(x) {
    mle <- 1/mean(x)
    
    return(list(MLE=mle, ll=sum(dstick(x, mle, log=TRUE)), df=1, nobs=length(x)))
}

.fitTnegb <- function(x) {
    fun <- function(par) {
        par[2] <- exp(par[2])
        return(-sum(dtnegb(x, par[1], par[2], log=TRUE)))
    }
    
    init.mu <- mean(x)
    if(var(x) < mean(x)) {
        init.k <- log(10000)
    } else {
        init.k <- log(mean(x)^2 / (var(x) - mean(x)))
    }
    
    fit <- optim(c(init.mu/2, init.k/2), fun, method='L-BFGS-B', lower=rep(.Machine$double.eps, 2), upper=c(init.mu*10, init.k*2))
    fit$par[2] <- exp(fit$par[2])
    
    # return(list(MLE=fit$par, ll=-fit$value, df=2, nobs=length(x)))
    return(fit)
}

.fitTpois <- function(x) {
    mle <- mean(x) + gsl::lambert_W0(-mean(x)*exp(-mean(x)))
    
    return(list(MLE=mle, ll=sum(dtpois(x, mle, log=TRUE)), df=1, nobs=length(x)))
}

