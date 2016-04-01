#' @title Log likelihood method for objects of class \code{sad}
#'  
#' @description Compute log likelihood for objects of class \code{sad}
#' 
#' @details See examples
#' 
#' @param x an object of class \code{sad}
#' 
#' @keywords species abundance, SAD
#' @export
#' 
#' @examples
#' x <- rtpois(1000, 2)
#' sad1 <- fitSAD(x, models='tpois')
#' logLik(sad1[[1]])
#' 
#' @return An object of class \code{logLik} with attributes "\code{df}" and "\code{nobs}".
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD
# @references

logLik.sad <- function(x) {
    out <- x$ll
    attr(out, 'df') <- x$df
    attr(out, 'nobs') <- x$nobs
    class(out) <- 'logLik'
    
    return(out)
}