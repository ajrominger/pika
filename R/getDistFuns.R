#' @title Functions to get distribution functions from objects of class \code{sad}
#'  
#' @description Get distribution functions from objects of class \code{sad}
#' 
#' @details See examples
#' 
#' @param s an object of class \code{sad}
#' @param ... arguments to be passed to distribution functions, e.g. \code{log}, \code{lower.tail}
#' 
#' @keywords species abundance, SAD
#' @export
#' 
#' @examples
#' x <- rtpois(1000, 2)
#' sad1 <- fitSAD(x, models='tpois')
#' getdfun(sad1[[1]])(1:10, log=TRUE)
#' 
#' @return A function (d, p, q, r) parameterized with the distribution of the \code{sad} object
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD
# @references

#' @export
#' @rdname getDistFuns
getdfun <- function(s) {
    switch(s$model, 
           'fish' = function(x, ...) dfish(x, s$MLE, ...),
           'plnorm' = function(x, ...) dplnorm(x, s$MLE[1], s$MLE[2], ...),
           'stick' = function(x, ...) dstick(x, s$MLE, ...),
           'tnegb' = function(x, ...) dtnegb(x, s$MLE[1], s$MLE[2], ...),
           'tpois' = function(x, ...) dtpois(x, s$MLE, ...))
}

#' @export
#' @rdname getDistFuns
getpfun <- function(s) {
    switch(s$model, 
           'fish' = function(q, ...) pfish(q, s$MLE, ...),
           'plnorm' = function(q, ...) pplnorm(q, s$MLE[1], s$MLE[2], ...),
           'stick' = function(q, ...) pstick(q, s$MLE, ...),
           'tnegb' = function(q, ...) ptnegb(q, s$MLE[1], s$MLE[2], ...),
           'tpois' = function(q, ...) ptpois(q, s$MLE, ...))
}


#' @export
#' @rdname getDistFuns
getqfun <- function(s) {
    switch(s$model, 
           'fish' = function(p, ...) qfish(p, s$MLE, ...),
           'plnorm' = function(p, ...) qplnorm(p, s$MLE[1], s$MLE[2], ...),
           'stick' = function(p, ...) qstick(p, s$MLE, ...),
           'tnegb' = function(p, ...) qtnegb(p, s$MLE[1], s$MLE[2], ...),
           'tpois' = function(p, ...) qtpois(p, s$MLE, ...))
}

#' @export
#' @rdname getDistFuns
getrfun <- function(s) {
    switch(s$model, 
           'fish' = function(n) rfish(n, s$MLE),
           'plnorm' = function(n) rplnorm(n, s$MLE[1], s$MLE[2]),
           'stick' = function(n) rstick(n, s$MLE),
           'tnegb' = function(n) rtnegb(n, s$MLE[1], s$MLE[2]),
           'tpois' = function(n) rtpois(n, s$MLE))
}
