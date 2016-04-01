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

