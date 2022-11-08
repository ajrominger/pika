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


## ================================
## function to simulate from an SAD
## ================================
#' @export

.simSAD <- function(x) {
    n <- length(x$data)
    r <- switch(x$model, 
                'lseries' = rlseries(n, x$MLE),
                'plnorm' = rplnorm(n, x$MLE[1], x$MLE[2]),
                'stick' = rstick(n, x$MLE),
                'tnegb' = rtnegb(n, x$MLE[1], x$MLE[2]),
                'tpois' = rtpois(n, x$MLE))
    
    return(r)
}
