context('distribution functions work')

test_that('truncated Poisson pfun works', {
    expect_true(all(ptpois(1:10, 1) - cumsum(dtpois(1:10, 1)) < .Machine$double.eps))
    expect_true(all(ptpois(1:10, 1, log=TRUE) - log(cumsum(dtpois(1:10, 1))) < .Machine$double.eps))
    expect_true(all(ptpois(1:10, 1, lower.tail=FALSE) - (1 - cumsum(dtpois(1:10, 1))) < .Machine$double.eps))
})
