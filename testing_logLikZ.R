library(socorro)
library(pika)

## ===================================================================
## script to see if we can accurately get the mean and sd for the null 
## distribution and if CLT holds up for small sample size
## ===================================================================

dfun <- function(x, ...) dfish(x, 0.1, ...)
S <- 100

foo <- replicate(1000, {
    x <- rfish(S, 0.1)
    sum(dfun(x, log = TRUE))
})

x0 <- 1:10^5
p0 <- dfun(1:10^5, log = TRUE)
p0 <- p0[is.finite(p0)]

m <- sum(p0 * exp(p0))
v <- sum((m - p0)^2 * exp(p0))

plot(simpECDF(foo))
curve(pnorm(x, m * S, sqrt(v * S)), col = 'red', add = TRUE)

mean(foo)
S * m

var(foo)
S * v


## =====================================================================
## script to explicitly test the z-value functions (old with simulation, 
## and new simulation free)
## =====================================================================

logLikZ_old <- function(x, nrep=1000, return.sim=FALSE) {
    lik.obs <- logLik(x)
    n <- x$nobs
    rfun <- getrfun(x)
    dfun <- getdfun(x)
    
    lik.sim <- replicate(nrep, {
        newx <- rfun(n)
        sum(dfun(newx, log=TRUE))
    })
    
    z <- ((lik.obs - mean(lik.sim)) / sd(lik.sim))^2
    
    if(return.sim) {
        lik.sim <- ((lik.sim - mean(lik.sim)) / sd(lik.sim))^2
    } else {
        lik.sim <- NULL
    }
    
    return(list(z=as.numeric(z), obs=lik.obs, sim=lik.sim))
}

logLikZ_new <- function(x) {
    lik.obs <- logLik(x)
    n <- x$nobs
    
    ## hypothetical distribution of probabilities
    p0 <- .p0(x)
    
    ## hypothetical mean and var
    m <- sum(p0 * exp(p0)) * n
    v <- sum((m/n - p0)^2 * exp(p0)) * n
    
    ## z^2-value
    z <- ((lik.obs - m) / sqrt(v))^2
    
    return(list(z=as.numeric(z), obs=lik.obs))
}

.p0 <- function(x) {
    n0 <- 1:10^5
    dfun <- getdfun(x)
    
    p0 <- dfun(n0, log = TRUE)
    p0 <- p0[is.finite(p0)]
    
    if(exp(p0[length(p0)]) > .Machine$double.eps^0.75) {
        n0add <- (n0[length(p0)] + 1):10^6
        p0add <- dfun(n0add, log = TRUE)
        p0add <- p0add[is.finite(p0add)]
        
        p0 <- c(p0, p0add)
    }
    
    return(p0)
}

x <- sad(rfish(100, 0.01), 'fish', keepData = TRUE)
plot(x, ptype = 'rad')
foo1 <- logLikZ_old(x)
foo2 <- logLikZ_new(x)

foo1$z
foo2$z
