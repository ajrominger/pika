#' @title Sample from an SAD object
#'  
#' @description Samples a desired number of individuals from an \code{sad} object.
#' 
#' @details If a single scalar values is supplied for \code{prob}, size will be 
#' calculated as \code{size = prob*sum(x)}; if a vector of length equal to length of \code{x} is giving this will
#' not change the value of \code{size} but rather indicate the relative probabilities of each species
#' appearing in the sample (e.g. for simulating species-specific detection bias)
#' 
#' @param x an object of class \code{sad}
#' @param size the number of individuals to sample from the SAD
#' @param replace logical, should sampling be done with replacement or without
#' @param prob an optional numeric specifying probability of inclusion in a sample, either a 
#' single scalar value or vector of length equal to length of \code{x}, see Details.
#' 
#' @keywords subsampling, species abundance, SAD
#' @export
#' 
#' @examples
#' x <- rtpois(1000, 2)
#' sad1 <- fitSAD(x, models='tpois')
#' sample(sad1[[1]], 10)
#' 
#' @return A numeric vector of variable length depending on how many species have sampled abundance > 0
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#' @seealso sample
# @references

sample.sad <- function(x, size, replace=FALSE, prob=NULL) {
    if(length(prob) >= 1) {
        prob[prob > 1] <- 1
        if(length(prob) == 1) size <- sum(x)*prob
    } else {
        prob <- 1
    }
    
    spp <- rep(1:length(x$data), x$data)
    
    samp <- sample(spp, size=size, replace=replace, prob=rep(prob, x$data))
    
    return(as.numeric(table(samp)))
}
