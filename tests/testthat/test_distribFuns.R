context('distribution functions work')

## =================================
## Truncated Poisson
## =================================

test_that('truncated Poisson dfun works', {
    expect_equal(dtpois(1:10, 1), dpois(1:10, 1) / (1 - dpois(0, 1)))
    expect_warning(dtpois(0.5, 1))
    expect_equal(dtpois(0, 1), 0)
})

test_that('truncated Poisson pfun works', {
    expect_true(all(abs(ptpois(1:10, 1) - cumsum(dtpois(1:10, 1))) < .Machine$double.eps))
    expect_true(all(abs(ptpois(1:10, 1, log=TRUE) - log(cumsum(dtpois(1:10, 1)))) < .Machine$double.eps))
    expect_true(all(abs(ptpois(1:10, 1, lower.tail=FALSE) - (1 - cumsum(dtpois(1:10, 1)))) < .Machine$double.eps))
    expect_equal(ptpois(0, 1), 0)
})

test_that('truncated Poisson qfun works', {
    expect_equal(qtpois(ptpois(1:30, 10), 10), 1:30)
    expect_equal(qtpois(ptpois(1:30, 10, lower.tail=FALSE), 10, lower.tail=FALSE), 1:30)
    expect_equal(qtpois(ptpois(1:30, 10, log=TRUE), 10, log=TRUE), 1:30)
})

test_that('truncated Poisson rfun works', {
    m <- 2
    r1 <- rtpois(10000, m)
    expect_true(abs(mean(r1) - m / (1 - dpois(0, m))) < 0.05)
})


## =================================
## Truncated Negative Binomial
## =================================

test_that('truncated Negative Binomial dfun works', {
    expect_equal(dtnegb(1:10, 2, 3), dnbinom(1:10, mu=2, size=3) / (1 - dnbinom(0, mu=2, size=3)))
    expect_warning(dtnegb(0.5, 2, 3))
    expect_equal(dtnegb(0, 2, 3), 0)
})

test_that('truncated Negative Binomial pfun works', {
    expect_true(all(abs(ptnegb(1:10, 2, 3) - cumsum(dtnegb(1:10, 2, 3))) <= .Machine$double.eps))
    expect_true(all(abs(ptnegb(1:10, 2, 3, log=TRUE) - log(cumsum(dtnegb(1:10, 2, 3)))) <= .Machine$double.eps*4))
    expect_true(all(abs(ptnegb(1:10, 2, 3, lower.tail=FALSE) - (1 - cumsum(dtnegb(1:10, 2, 3)))) <= .Machine$double.eps))
    expect_equal(ptnegb(0, 2, 3), 0)
})

test_that('truncated Negative Binomial qfun works', {
    expect_equal(qtnegb(ptnegb(1:10, 2, 3), 2, 3), 1:10)
    expect_equal(qtnegb(ptnegb(1:10, 2, 3, lower.tail=FALSE), 2, 3, lower.tail=FALSE), 1:10)
    expect_equal(qtnegb(ptnegb(1:10, 2, 3, log=TRUE), 2, 3, log=TRUE), 1:10)
})

test_that('truncated Negative Binomial rfun works', {
    m <- 2
    r1 <- rtnegb(10000, m, 3)
    expect_true(abs(mean(r1) - m / (1 - dnbinom(0, mu=m, size=3))) < 0.05)
})

## =================================
## Fisher log series
## =================================

test_that('Fisher log series dfun works', {
    expect_true(all(abs(VGAM::dlog(1:10, exp(-0.1)) - dfish(1:10, 0.1)) <= .Machine$double.eps))
    expect_warning(dfish(0.5, 0.1))
    expect_equal(dfish(0, 0.1), 0)
})

test_that('Fisher log series pfun works', {
    expect_true(all(abs(pfish(1:10, 0.1) - cumsum(dfish(1:10, 0.1))) <= .Machine$double.eps*6))
    expect_true(all(abs(pfish(1:10, 0.1, log=TRUE) - log(cumsum(dfish(1:10, 0.1)))) <= .Machine$double.eps*12))
    expect_true(all(abs(pfish(1:10, 0.1, lower.tail=FALSE) - (1 - cumsum(dfish(1:10, 0.1)))) <= .Machine$double.eps*6))
    expect_equal(pfish(0, 0.1), 0)
})

test_that('Fisher log series qfun works', {
    expect_equal(qfish(pfish(1:10, 0.1), 0.1), 1:10)
    expect_equal(qfish(pfish(1:10, 0.1, lower.tail=FALSE), 0.1, lower.tail=FALSE), 1:10)
    expect_equal(qfish(pfish(1:10, 0.1, log=TRUE), 0.1, log=TRUE), 1:10)
})

test_that('Fisher log series rfun works', {
    b <- 0.1
    tru.mean <- -1/log(1-exp(-b)) * exp(-b)/(1 - exp(-b))
    r1 <- rfish(10000, b)
    tru.mean - mean(r1)
    expect_true(abs(mean(r1) - tru.mean) < 0.5)
})

## =================================
## Broken stick
## =================================

test_that('Broken stick dfun works', {
    expect_equal(dstick(1:10, 0.5), dgeom(1:10, 0.5) / (1 - dgeom(0, 0.5)))
    expect_warning(dstick(0.5, 1))
    expect_equal(dstick(0, 1), 0)
})

test_that('Broken stick pfun works', {
    expect_true(all(abs(pstick(1:10, 0.5) - cumsum(dstick(1:10, 0.5))) < .Machine$double.eps))
    expect_true(all(abs(pstick(1:10, 0.5, log=TRUE) - log(cumsum(dstick(1:10, 0.5)))) < .Machine$double.eps))
    expect_true(all(abs(pstick(1:10, 0.5, lower.tail=FALSE) - (1 - cumsum(dstick(1:10, 0.5)))) < .Machine$double.eps))
    expect_equal(pstick(0, 1), 0)
})

test_that('Broken stick qfun works', {
    expect_equal(qstick(pstick(1:10, 0.5), 0.5), 1:10)
    expect_equal(qstick(pstick(1:10, 0.5, lower.tail=FALSE), 0.5, lower.tail=FALSE), 1:10)
    expect_equal(qstick(pstick(1:10, 0.5, log=TRUE), 0.5, log=TRUE), 1:10)
})

test_that('Broken stick rfun works', {
    p <- 0.5
    r1 <- rstick(10000, p)
    expect_true(abs(mean(r1) - 1/p) < 0.05)
})
