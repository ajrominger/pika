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

.fitFish <- function(x) {
    xbar <- mean(x)
    # fun <- function(p) -1/log(1 - p) * p/(1 - p) - xbar
    # uniroot(fun, lower=.Machine$double.eps, upper=1-.Machine$double.eps, tol=.Machine$double.eps)
    fun <- function(b) log((-1/log(1 - exp(-b)) * exp(-b)/(1 - exp(-b)) - xbar)^2)
    fit <- optimize(fun, lower=xbar^-4, upper=20/xbar, tol=.Machine$double.eps)
    
    ## consider doing log(b)
}


# x <- rfish(100, 0.01)
.fitFish(x)

