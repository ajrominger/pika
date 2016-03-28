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
