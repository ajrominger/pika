context('fitting functions work')

## =================================
## General functionality
## =================================

test_that('fitting function returns expected results', {
    fitAll <- fitSAD(1:10)
    fitOne <- fitSAD(1:10, models='fish')
    expect_equal(length(fitAll), 5)
    expect_is(fitOne[[1]], 'sad')
    expect_is(fitOne, 'list')
    expect_equal(names(fitOne), 'fish')
})


## =================================
## Fisher
## =================================

test_that('Fisher log series fitting works', {
    b <- 0.01
    x <- rfish(10000, b)
    fit <- fitSAD(x, 'fish')
    expect_true(abs(fit[[1]]$MLE - b) < 0.05*b)
})


## =================================
## Poisson Log Norm
## =================================

test_that('Poisson log normal fitting works', {
    m <- 1
    s <- 3
    x <- rplnorm(10000, m, s)
    fit <- fitSAD(x, 'plnorm')
    expect_true(abs(fit[[1]]$MLE[1] - m) < 0.05*m)
    expect_true(abs(fit[[1]]$MLE[2] - s) < 0.05*s)
})


## =================================
## Broken Stick
## =================================

test_that('Broken stick fitting works', {
    r <- 0.8
    x <- rstick(10000, r)
    fit <- fitSAD(x, 'stick')
    expect_true(abs(fit[[1]]$MLE - r) < 0.05*r)
})


## =================================
## Truncated Neg Binomial
## =================================

test_that('Truncated negative binomial fitting works', {
    m <- 1
    k <- 3
    x <- rtnegb(10000, m, k)
    fit <- fitSAD(x, 'tnegb')
    expect_true(abs(fit[[1]]$MLE[1] - m) < 0.05*m)
    expect_true(abs(fit[[1]]$MLE[2] - k) < 0.05*k)
})


## =================================
## Truncated Poisson
## =================================

test_that('truncated Poisson fitting works', {
    la <- 3
    x <- rtpois(10000, la)
    fit <- fitSAD(x, 'tpois')
    expect_true(abs(fit[[1]]$MLE - la) < 0.05*la)
})
