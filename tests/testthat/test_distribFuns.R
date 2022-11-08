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
## log series
## =================================

test_that('log series dfun works', {
    expect_true(all(abs(VGAM::dlog(1:10, exp(-0.1)) - dlseries(1:10, 0.1)) <= .Machine$double.eps))
    expect_warning(dlseries(0.5, 0.1))
    expect_equal(dlseries(0, 0.1), 0)
})

test_that('log series pfun works', {
    expect_true(all(abs(plseries(1:10, 0.1) - cumsum(dlseries(1:10, 0.1))) <= .Machine$double.eps*6))
    expect_true(all(abs(plseries(1:10, 0.1, log=TRUE) - log(cumsum(dlseries(1:10, 0.1)))) <= .Machine$double.eps*12))
    expect_true(all(abs(plseries(1:10, 0.1, lower.tail=FALSE) - (1 - cumsum(dlseries(1:10, 0.1)))) <= .Machine$double.eps*6))
    expect_equal(plseries(0, 0.1), 0)
})

test_that('log series qfun works', {
    expect_equal(qlseries(plseries(1:10, 0.1), 0.1), 1:10)
    expect_equal(qlseries(plseries(1:10, 0.1, lower.tail=FALSE), 0.1, lower.tail=FALSE), 1:10)
    expect_equal(qlseries(plseries(1:10, 0.1, log=TRUE), 0.1, log=TRUE), 1:10)
})

test_that('log series rfun works', {
    b <- 0.1
    tru.mean <- -1/log(1-exp(-b)) * exp(-b)/(1 - exp(-b))
    r1 <- rlseries(10000, b)
    tru.mean - mean(r1)
    expect_true(abs(mean(r1) - tru.mean) < 0.5)
})

## =================================
## Poisson log normal
## =================================

test_that('Poisson log normal dfun works', {
    expect_equal(dplnorm(1:10, 1, 1), poilog::dpoilog(1:10, 1, 1) / (1 - poilog::dpoilog(0, 1, 1)))
    expect_warning(dplnorm(0.5, 1, 1))
    expect_equal(dplnorm(0, 1, 1), 0)
})

test_that('Poisson log normal pfun works', {
    expect_true(all(abs(pplnorm(1:10, 1, 1) - cumsum(dplnorm(1:10, 1, 1))) < .Machine$double.eps))
    expect_true(all(abs(pplnorm(1:10, 1, 1, log=TRUE) - log(cumsum(dplnorm(1:10, 1, 1)))) < .Machine$double.eps))
    expect_true(all(abs(pplnorm(1:10, 1, 1, lower.tail=FALSE) - (1 - cumsum(dplnorm(1:10, 1, 1)))) < .Machine$double.eps))
    expect_equal(pplnorm(0, 1, 1), 0)
})

test_that('Poisson log normal qfun works', {
    expect_equal(qplnorm(pplnorm(1:10, 1, 1), 1, 1), 1:10)
#    expect_equal(qplnorm(pplnorm(1:10, 1, 1, lower.tail=FALSE), 1, 1, lower.tail=FALSE), 1:10)
    expect_equal(qplnorm(pplnorm(1:10, 1, 1, log=TRUE), 1, 1, log=TRUE), 1:10)
})

test_that('Poisson log normal rfun works', {
    m <- 2
    s <- 1
    r1 <- rplnorm(10000, m, s)
    expect_true(abs(mean(r1) - exp(m+s^2/2) / (1 - poilog::dpoilog(0, m, s))) < 0.5)
})
