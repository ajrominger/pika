#' @title Fisher's logseries
#'  
#' @description \code{dfish} gives the probability mass function, \code{pfish} gives the cumulative mass 
#' function, \code{qfish} the quantile function, \code{rfish} randome number generation
#' 
#' @details These functions assume infinite support of the Fisher logseries
#' 
#' @param x
#' @param beta
#' @param log
#' @param n
#' 
#' @keywords Fisher logseries, species abundance, SAD
#' @export
#' 
#' @examples
#' 
#' dfish(1:10, 0.01)
#' 
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso 
# @references 

dfish <- function(x, beta, log=FALSE) {
    out <- 1/log(1/(1-exp(-beta))) * exp(-beta*x)/x
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer x = %s', bad))
        }
        
        out[x %% 1 != 0] <- 0
    }
    
    if(log) out <- log(out)
    return(out)
}


#' @rdname dfish

pfish <- function(q, beta, lower.tail=TRUE, log=FALSE) {
    out <- 1 + .betax(exp(-beta), q+1, 0) / log(1 - exp(-beta))
    
    if(any(q %% 1 != 0)) {
        for(bad in q[q %% 1 != 0]) {
            warning(sprintf('non-integer q = %s', bad))
        }
        
        out[q %% 1 != 0] <- 0
    }
    
    if(lower.tail) out <- 1 - out
    if(log) out <- log(out)
    return(out)
}

#' @rdname dfish

qfish <- function(p, beta, lower.tail=TRUE, log=FALSE) {
    approx(x=.fishcdf(1:10000, beta), y=1:10000, xout=p, method='constant', yleft=NaN, yright=NaN, f=0)
}

## =================================
## helper functions
## =================================

.betax <- function(x, a, b) {
    1/a * x^a * (1-x)^b * gsl::hyperg_2F1(a+b, 1, a+1, x)
}

.fishcdf <- function(x, beta) {
    1 + .betax(exp(-beta), x+1, 0) / log(1 - exp(-beta))
}

