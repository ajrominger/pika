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
                'fish' = rfish(n, x$MLE, log=log),
                'plnorm' = rplnorm(n, x$MLE[1], x$MLE[2], log=log),
                'stick' = rstick(n, x$MLE, log=log),
                'tnegb' = rtnegb(n, x$MLE[1], x$MLE[2], log=log),
                'tpois' = rtpois(n, x$MLE, log=log))
    
    return(r)
}
