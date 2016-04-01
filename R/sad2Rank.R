
#' @title meteDist2Rank
#'  
#' @description Calculate the rank distribution of a sad object 
#'
#' @details
#' Extracts the predicted rank distribution from a \code{meteDist} object. 
#' This is effectively the quantile function of the distribution. Used, e.g., 
#' in \code{plot.meteDist}
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
#' @seealso fitSAD
#' @references 

sad2Rank <- function(x, S=NULL) {
    if(is.null(S)) {
        S <- length(x$data)
    }
    
    qfun <- switch(x$model,
                   'fish' = function(p) qfish(p, x$MLE),
                   'plnorm' = function(p) qplnorm(p, x$MLE[1], x$MLE[2]),
                   'stick' = function(p) qstick(p, x$MLE),
                   'tnegb' = function(p) qtnegb(p, x$MLE[1], x$MLE[2]),
                   'tpois' = function(p) qtpois(p, x$MLE))
    
    return(qfun(seq(1, 1/S, length=S) - 1/(2*S)))
}

