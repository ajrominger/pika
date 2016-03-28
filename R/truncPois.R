#' @title Truncated Poisson distribuiton
#'  
#' @description \code{dtpois} gives the probability mass function, \code{ptpois} gives the cumulative mass 
#' function, \code{qtpois} the quantile function, \code{rtpois} randome number generation
#' 
#' @details These functions assume infinite support of the Poisson log-normal from [1, Inf).
#' 
#' @param x vector of integers for which to return the probability
#' @param q vector of integers for which to return the cumulative probability
#' @param p vector of probabilities for which to return the quantile
#' @param n number of random replicates
#' @param lambda mean abundance, see Details
#' @param log logical, should the log probability be used
#' @param lower.tail logical, should the lower tail be used
#' 
#' @keywords Truncated Poisson, species abundance, SAD
#' @export
#' 
#' @examples
#' 
#' dtpois(1:10, 2)
#' 
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso 
# @references 

#' @rdname TPois

dtpois <- function(x, lambda, log=FALSE) {
    out <- dpois(x, lambda) / (1 - dpois(0, lambda))
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer x = %s', bad))
        }
        
        out[x %% 1 != 0] <- 0
    }
    
    if(log) out <- log(out)
    return(out)
}


#' @rdname TPois

ptpois <- function(q, lambda, lower.tail=TRUE, log=FALSE) {
    if(log) {
#         out <- ppois(q, lambda, lower.tail=lower.tail, log=TRUE) - (1 - dpois(0, lambda, log=TRUE))
        out <- NA
    } else {
        out <- (ppois(q, lambda, lower.tail=lower.tail) - dpois(0, lambda)) / (1 - dpois(0, lambda))
    }
    
    return(out)
}

#' @rdname TPois

qtpois <- function(p, lambda, lower.tail=TRUE, log=FALSE) {
#     if(log) p <- exp(p)
#     if(!lower.tail) p <- 1 - p
#     
#     out <- approx(x=cumsum(dnbinom(1:10000, mu=mu, size=k)) / (1 - dnbinom(0, mu=mu, size=k)), y=1:10000,
#                   xout=p, xout=p, method='constant', yleft=NaN, yright=NaN, f=0)
#     
#     if(any(is.nan(out))) {
#         warning('NaNs produced')
#     }
#     
#     return(out)
    qpois(p/(1-dpois(0, lambda)), lambda)
}


#' @rdname TPois

rtpois <- function(n, lambda) {
    N <- 100 * n / (1 - dpois(0, lambda))
    temp <- rpois(N, lambda)
    temp <- temp[temp > 0]
    
    if(length(temp) < n) warning(sprintf('could not find %s unique random variates, using bootstrapping', n))
    
    return(sample(temp, n, rep=ifelse(n < length(temp), FALSE, TRUE)))
}
