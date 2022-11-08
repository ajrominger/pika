#' @title Fit a model for a species abundance distribution
#'  
#' @description \code{fitSAD} uses maximum likelihood to fit one or several SAD models to data
#' 
#' @details The default behavior fits all models availible in \code{pika}:
#' \describe{
#'  \item{\code{fish}}{Fisher's log series}
#'  \item{\code{plnorm}}{Poisson log normal}
#'  \item{\code{stick}}{Broken stick}
#'  \item{\code{tnegb}}{Truncated negative binomial}
#'  \item{\code{tpois}}{Truncated Poisson}
#' }
#' Exact solutions are availible for \code{tpois} and \code{stick}. In the case of \code{fish}, 
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
#' @seealso logLik.sad, dfish, dplnorm, dstick, dtnegb, dtpois
# @references 

fitSAD <- function(x, models=c('fish', 'plnorm', 'stick', 'tnegb', 'tpois'), keepData=FALSE) {
    if(length(models) > 1) {
        out <- lapply(models, function(m) fitSAD(x, m, keepData)[[1]])
        # return(out)
    } else {
        out <- switch(models,
                      'fish' = .fitFish(x),
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

