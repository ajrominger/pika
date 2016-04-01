#' @title Print method for objects of class \code{sad}
#'  
#' @description Print objects of class \code{sad}
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
#' print(sad1[[1]])
#' 
#' @return Silently returns \code{x}.
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD
# @references

print.sad <- function(x) {
    cat(sprintf('species abundance distribution modeled by "%s" with parameter%s \n', x$model, ifelse(length(x$MLE > 1), '', 's')))
    print(x$MLE)
    cat(sprintf('%s data', ifelse(is.null(x$data), 'does not include', 'includes')))
    invisible(x)
}



#' @title Plot method for objects of class \code{sad}
#'  
#' @description Plot objects of class \code{sad}
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
#' plot(sad1[[1]])
#' 
#' @return \code{NULL}.
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso fitSAD
# @references

plot.sad <- function(x, ptype=c('cdf', 'rad'), ...) {
    ptype <- match.arg(ptype, c('cdf', 'rad'))
    
    if(!is.null(x$data)) {
        if(ptype=='cdf') {
            thisCDF <- .ecdf(x$data)
            plot(thisCDF[, 1], thisCDF[, 2], ...)
        } else {
            plot(sort(x$data, decreasing=TRUE))
        }
    }
    
#     if(!is.null(x$model)) {
#         if(dev.cur() == 1) {
#             plotfun <- function(x, y) plot(x, y, ...)
#         } else {
#             plotfun <- function(x, y) lines(x, y)
#         }
#         if(ptype)
#     }
}


## =================================
## simplified empirical CDF function
## =================================
#' @export
.ecdf <- function(x) {
    tabz <- table(x)
    yval <- cumsum(as.numeric(tabz))/sum(as.numeric(tabz))
    xval <- as.numeric(names(tabz))
    
    return(cbind(xval, yval))
}