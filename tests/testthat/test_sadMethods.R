context('sad methods work')

## =================================
## sad2Rank
## =================================

test_that('sad2Rank works', {
    la <- 2
    x <- rtpois(100, la)
    xfit <- fitSAD(x, 'tpois', keepData=TRUE)[[1]]
    r <- sad2Rank(xfit)
    expect_equal(length(r), length(x))
    expect_true(abs(fitSAD(r, 'tpois')[[1]]$MLE - la) < 0.05*la)
})


## =================================
## logLik
## =================================

test_that('logLik works', {
    la <- 2
    x <- rtpois(100, la)
    xfit <- fitSAD(x, 'tpois', keepData=TRUE)[[1]]
    expect_is(logLik(xfit), 'logLik')
    z <- logLikZ(xfit, nrep=10, return.sim=TRUE)
    expect_equal(length(z), 3)
    expect_equal(z$obs, logLik(xfit))
    expect_equal(length(z$sim), 10)
})


## =================================
## residuals and mse
## =================================

test_that('residuals and mse work', {
    la <- 2
    x <- rtpois(100, la)
    xfit <- fitSAD(x, 'tpois', keepData=TRUE)[[1]]
    resid1 <- residuals(xfit)
    resid2 <- residuals(xfit, type='cumulative')
    resid3 <- residuals(xfit, log=TRUE, relative=FALSE)
    resid4 <- residuals(xfit, log=FALSE, relative=TRUE)
    resid5 <- residuals(xfit, type='cumulative', log=TRUE, relative=TRUE)
    resid6 <- residuals(xfit, type='cumulative', log=FALSE, relative=TRUE)
    
    expect_equal(length(resid1), 100)
    expect_equal(mse(xfit), mean(resid1^2))
    expect_equal(mse(xfit, type='cumulative'), mean(resid2^2))
    expect_equal(mse(xfit, log=TRUE, relative=FALSE), mean(resid3^2))
    expect_equal(mse(xfit, log=FALSE, relative=TRUE), mean(resid4^2))
    expect_equal(mse(xfit, type='cumulative', log=TRUE, relative=TRUE), mean(resid5^2))
    expect_equal(mse(xfit, type='cumulative', log=FALSE, relative=TRUE), mean(resid6^2))
    
    z <- mseZ(xfit, nrep=10, return.sim=TRUE, type='cumulative', log=FALSE, relative=TRUE)
    expect_equal(length(z), 3)
    expect_equal(z$obs, mse(xfit, type='cumulative', log=FALSE, relative=TRUE))
    expect_equal(length(z$sim), 10)
})



## =================================
## get distrib funs
## =================================

test_that('get distrib funs works', {
    x <- rtpois(100, 2)
    xfit <- fitSAD(x, 'tpois', keepData=TRUE)[[1]]
    expect_is(getdfun(xfit), 'function')
    expect_equal(getdfun(xfit)(1:10), dtpois(1:10, xfit$MLE))
    expect_equal(getpfun(xfit)(1:10), ptpois(1:10, xfit$MLE))
    expect_equal(getqfun(xfit)(1/(1:10)), qtpois(1/(1:10), xfit$MLE))
    expect_equal(length(getrfun(xfit)(10)), 10)
})
