#' @title Fisher's logseries
#'  
#' @description \code{dfish} gives the probability mass function, \code{pfish} gives the cumulative mass 
#' function, \code{qfish} the quantile function, \code{rfish} randome number generation
#' 
#' @details
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
#' @seealso 
#' @references 

dfish <- function(x, beta, log=FALSE) {
    out <- 1/log(1/(1-exp(-beta))) * exp(-beta*x)/x
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer n = %s', bad))
        }
        
        out[x %% 1 != 0] <- 0
    }
    
    if(log) out <- log(out)
    return(out)
}

pfish <- function(n, beta, log=FALSE) {
    out <- 1 + .betax(exp(-beta), n+1, 0) / log(1 - exp(-beta))
    
    if(any(n %% 1 != 0)) {
        for(bad in n[n %% 1 != 0]) {
            warning(sprintf('non-integer n = %s', bad))
        }
        
        out[n %% 1 != 0] <- 0
    }
    
    if(log) out <- log(out)
    return(out)
}

.betax <- function(x, a, b) {
    1/a * x^a * (1-x)^b * gsl::hyperg_2F1(a+b, 1, a+1, x)
}

.fishcdf <- function(n, beta) {
    1 + .betax(exp(-beta), n+1, 0) / log(1 - exp(-beta))
}
