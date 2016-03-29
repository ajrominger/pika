## ===============================================
## function implementing inverse tranform sampling
## for descrete random number generation
## ===============================================

.its <- function(n, qfun, minp, ...) {
    r <- runif(n, minp, 1)
#     browser()
    return(floor(qfun(r, ...)))
}
