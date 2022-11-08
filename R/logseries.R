#' @title Logseries distribution
#'  
#' @description \code{dlseries} gives the probability mass function, \code{plseries} gives the cumulative mass 
#' function, \code{qlseries} the quantile function, \code{rlseries} randome number generation
#' 
#' @details These functions assume infinite support of the logseries from [1, Inf).  An alternate parameterization is used to allow a larger range (i.e. the positive number line) of the parameter \code{beta}: \eqn{\frac{1}{log(1/(1-e^{-\beta}))} * \frac{e^{-\beta x}}{x}}.
#' 
#' @param x vector of integers for which to return the probability
#' @param q vector of integers for which to return the cumulative probability
#' @param p vector of probabilities for which to return the quantile
#' @param n number of random replicates
#' @param beta parameter of the log series, see Details
#' @param log logical, should the log probability be used
#' @param lower.tail logical, should the lower tail be used
#' 
#' @keywords logseries, species abundance, SAD
#' @export
#' 
#' @examples
#' 
#' dlseries(1:10, 0.01)
#' 
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
# @seealso 
# @references 

#' @rdname logseries

dlseries <- function(x, beta, log=FALSE) {
    out <- 1/log(1/(1-exp(-beta))) * exp(-beta*x)/x
    
    if(any(x %% 1 != 0)) {
        for(bad in x[x %% 1 != 0]) {
            warning(sprintf('non-integer x = %s', bad))
        }
        
        out[x %% 1 != 0] <- 0
    }
    
    out[x < 1] <- 0
    
    if(log) out <- log(out)
    return(out)
}


#' @export
#' @rdname logseries

plseries <- function(q, beta, lower.tail=TRUE, log=FALSE) {
    extraDistr::plgser(q, theta = exp(-beta), lower.tail = lower.tail, 
                       log.p = log)
}


#' @export
#' @rdname logseries

qlseries <- function(p, beta, lower.tail=TRUE, log=FALSE) {
    extraDistr::qlgser(p, theta = exp(-beta), lower.tail = lower.tail, 
                       log.p = log)
}


#' @export
#' @rdname logseries

rlseries <- function(n, beta) {
    extraDistr::rlgser(n, exp(-beta))
}

