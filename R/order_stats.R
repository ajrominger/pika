dorder <- function(x, r, s) {
    n <- s$nobs
    k <- n - r + 1
    dfun <- getdfun(s)
    pfun <- getpfun(s)
    
    jj <- 0:(n - k)
    
    nChooseJ <- choose(n, jj)
    p <- pfun(x)
    d <- dfun(x)
    
    a1 <- outer(1 - p, jj, '^')
    a2 <- outer(p, n - jj, '^')
    a3 <- outer(1 - p + d, jj, '^')
    a4 <- outer(p - d, n - jj, '^')
    
    o <- as.vector((a1 * a2 - a3 * a4) %*% nChooseJ)
    
    o[o < .Machine$double.eps^0.75] <- 0
    
    return(o)
}


porder <- function(q, r, s) {
    # create all objects that only depend on s
    n <- s$nobs
    pfun <- getpfun(s)
    
    # now objects that depend on s and q
    p <- pfun(q)
    
    # now determine if multiple of single r given
    if(length(r) == 1) {
        k <- n - r + 1
        jj <- 0:(n - k)
        nChooseJ <- choose(n, jj)
        
        a1 <- outer(1 - p, jj, '^')
        a2 <- outer(p, n - jj, '^')
        
        o <- as.vector((a1 * a2) %*% nChooseJ)
        
        o[o < .Machine$double.eps^0.75] <- 0
    } else {
        # multiple r given
        
        if(length(p) == 1) {
            p <- rep(p, length(r))
        } else if(length(p) != length(r)) {
            stop('lengths of `q` and `r` must be equal or one must be of length 1')
        }
        
        o <- sapply(1:length(r), function(i) {
            k <- n - r[i] + 1
            jj <- 0:(n - k)
            nChooseJ <- choose(n, jj)
            
            sum(nChooseJ * p[i]^jj * (1 - p[i])^(n - jj))
        })
        
        o[o < .Machine$double.eps^0.75] <- 0
    } 
    
    return(o)
}


.porder <- function(q, s) {
    q <- sort(q, decreasing = TRUE)
    # create all objects that only depend on s
    n <- s$nobs
    nChooseJ <- choose(n, 0:(n - 1))
    pfun <- getpfun(s)
    
    # now objects that depend on s and q
    p <- pfun(q)
    
    # now objects that depend on r
    o <- sapply(1:n, function(k) {
        ii <- 1:(n - k + 1)
        jj <- ii - 1
        
        1 - sum(nChooseJ[ii] * p[k]^(jj) * (1 - p[k])^(n - jj))
    })
    
    o[o < .Machine$double.eps^0.75] <- 0
    
    return(o)
}


Ps <- function(s) {
    n <- s$nobs
    ii <- 1:n
    pii <- .porder(s$data, s)
    
    o <- -n - 1/n * sum((2 * n + 1 - 2 * ii) * log(pii) +
                            (2 * ii - 1) * log(1 - pii))
    
    return(o)
}


Qs <- function(s) {
    pii <- .porder(s$data, s)
    
    return(sqrt(mean((pii - 0.5)^2)))
}




library(socorro)
library(viridis)

r <- rtnegb(500, 20, 2)
s <- sad(r, 'tnegb', keepData = TRUE)
plot(s, ptype = 'rad', log = 'y')

Qs(s)

plot(.porder(s$data, s))

sad2Rank(s)
nn <- 1:5000

ordRad <- lapply(1:s$nobs, function(i) {
    xs <- dorder(nn, i, s)
    Fs <- cumsum(xs)
    nmean <- sum(nn * xs)
    nmax <- which.max(xs)
    nmid <- min(which(Fs > 0.5))
    
    return(c(nmean = nmean, nmax = nmax, nmid = nmid))
})

ordRad <- do.call(rbind, ordRad)

plot(sad2Rank(s), log = 'y')
lines(ordRad[, 1], col = 'red')
lines(ordRad[, 2], col = 'blue')
lines(ordRad[, 3], col = 'purple')


plot(sad2Rank(s), ordRad[, 3], log = 'xy'); abline(0, 1)

Ps(s)

pp <- .porder(s$data, s)

layout(matrix(1:2, nrow = 1))
plot(sort(s$data, TRUE), col = quantCol(pp, viridis(20)), log = 'y')
lines(sad2Rank(s))

plot(sort(pp), col = quantCol(sort(pp), viridis(20)))

x <- sort(s$data, TRUE)
pp <- sapply(1:length(x), function(i) {
    sum(dorder(1:(x[i]), i, s))
})



plot(pp, .porder(s$data, s))
abline(0, 1)

porder(10:1, 81:90, s)
porder(1:10, 81, s)


Ps <- function(s) {
    n <- s$nobs
    ii <- 1:n
    
    pi <- sapply(ii, function(i) {
        2 * n + 1 - 2 * i * log(porder(x[i], i, s))
    })
}