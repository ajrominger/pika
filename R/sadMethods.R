

print.sad <- function(x) {
    cat(sprintf('species abundance distribution modeled by "%s" with parameters \n', x$model))
    print(x$MLE)
}

logLik.sad <- function(x) {
    out <- x$ll
    attr(out, 'df') <- x$df
    attr(out, 'nobs') <- x$nobs
    class(out) <- 'logLik'
    
    return(out)
}

plot.sad <- function(x, ptype=c('cdf', 'rad'), ...) {
    ptype <- match.arg(ptype, c('cdf', 'rad'))
    
    if(!is.null(x$data)) {
        if(ptype=='cdf') {
            tabz <- table(x$data)
            yval <- cumsum(as.numeric(tabz))/sum(as.numeric(tabz))
            xval <- as.numeric(names(tabz))
            
            plot(xval, yval, ...)
        } else {
            plot(sort(x$data, decreasing=TRUE))
        }
    } else {
        
    }
    
    if(!is.null(x$model)) {
        if(ptype)
    }
}