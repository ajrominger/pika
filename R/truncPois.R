#' @title Truncated Poisson distribuiton
#'  
#' @description \code{dtpois} gives the probability mass function, \code{ptpois} gives the cumulative mass 
#' function, \code{qtpois} the quantile function, \code{rtpois} randome number generation
#' 
#' @details These functions assume infinite support of the truncated Poisson from [1, Inf).
#' 
#' @param x vector of integers for which to return the probability
#' @param q vector of integers for which to return the cumulative probability
#' @param p vector of probabilities for which to return the quantile
#' @param n number of random replicates
#' @param lambda mean abundance
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
#' @seealso dpois
# @references 

#' @rdname TPois

dtpois <- function(x, lambda, log=FALSE) {
    if(log) {
        out <- dpois(x, lambda, log=TRUE) - log(1 - dpois(0, lambda))
        out[x < 1] <- -Inf
    } else {
        out <- dpois(x, lambda) / (1 - dpois(0, lambda))
        out[x < 1] <- 0
    }
    
    return(out)
}


#' @export
#' @rdname TPois

ptpois <- function(q, lambda, lower.tail=TRUE, log=FALSE) {
    if(lower.tail) {
        out <- (ppois(q, lambda, lower.tail=lower.tail) - dpois(0, lambda)) / (1 - dpois(0, lambda))
    } else {
        out <- ppois(q, lambda, lower.tail=lower.tail) / (1 - dpois(0, lambda))
    }
    
    out[q < 1] <- 0
    if(log) out <- log(out)
    
    return(out)
}


#' @export
#' @rdname TPois

qtpois <- function(p, lambda, lower.tail=TRUE, log=FALSE) {
    if(log) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    
    return(qpois(p*(1-dpois(0, lambda)) + dpois(0, lambda), lambda))
}


#' @export
#' @rdname TPois

rtpois <- function(n, lambda) {
    return(.its(n, qtpois, dtpois(0, lambda), lambda=lambda))
}
