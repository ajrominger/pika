context('fitting functions work')

## =================================
## Truncated Poisson
## =================================

test_that('truncated Poisson fitting works', {
    expect_equal(dtpois(1:10, 1), dpois(1:10, 1) / (1 - dpois(0, 1)))
    expect_warning(dtpois(0.5, 1))
    expect_equal(dtpois(0, 1), 0)
})