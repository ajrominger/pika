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

test_that('logLik.sad works', {
    la <- 2
    x <- rtpois(100, la)
    xfit <- fitSAD(x, 'tpois', keepData=TRUE)[[1]]
    expect_is(logLik(xfit), 'logLik')
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
