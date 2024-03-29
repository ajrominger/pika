
#' @title meteDist2Rank
#'  
#' @description Calculate the rank distribution of a sad object 
#'
#' @details
#' Extracts the predicted rank distribution from a \code{sad} object. 
#' This is effectively the quantile function of the distribution.
#' 
#' @param x \code{meteDist} object
# @keywords 
#' @export
#' 
#' @examples
#' 'stub'
#'                
#' @return A vector of predicted ranked abundances
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD, qqplot
# @references 

sad2Rank <- function(x, S=NULL) {
    if(is.null(S)) {
        S <- x$nobs
        if(S == 0 | is.null(S)) stop('nobs or data not provided in sad object "x"')
    }
    
    qfun <- switch(x$model,
                   'lseries' = function(p) qlseries(p, x$MLE),
                   'plnorm' = function(p) qplnorm(p, x$MLE[1], x$MLE[2]),
                   'stick' = function(p) qstick(p, x$MLE),
                   'tnegb' = function(p) qtnegb(p, x$MLE[1], x$MLE[2]),
                   'tpois' = function(p) qtpois(p, x$MLE))
    
    return(qfun(seq(1, 1/S, length=S) - 1/(2*S)))
}
