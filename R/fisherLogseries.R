## ==============================================
## distribution functions for Fisher's log series
## re-parameterized as in Harte and others
## ==============================================

dfish <- function(n, beta, log=FALSE) {
    out <- 1/log(1/(1-exp(-beta))) * exp(-beta*n)/n
    
    if(any(n %% 1 != 0)) {
        for(bad in n[n %% 1 != 0]) {
            warning(sprintf('non-integer n = %s', bad))
        }
        
        out[n %% 1 != 0] <- 0
    }
    
    if(log) out <- log(out)
    return(out)
}

pfish <- function(n, beta, log=FALSE)
