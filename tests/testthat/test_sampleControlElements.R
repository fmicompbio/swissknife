test_that("sampleControlElements() works properly", {
    i.control <- seq.int(10000)
    i.target <- seq.int(500) + length(i.control)
    x <- c(runif(length(i.control), min = 0, max = 10),
           rnorm(length(i.target), mean = 5, sd = 1))
    set.seed(1)
    expect_true(is.numeric(s1 <- sampleControlElements(x, idxTarget = i.target, idxControl = i.control)))
    set.seed(1)
    expect_true(is.numeric(s2 <- sampleControlElements(x,
                                      idxTarget = rep(c(FALSE,TRUE),c(length(i.control),length(i.target))),
                                      idxControl = NULL)))
    set.seed(1)
    expect_true(is.numeric(s3 <- sampleControlElements(list(x),
                                      idxTarget = rep(c(FALSE,TRUE),c(length(i.control),length(i.target))),
                                      idxControl = i.control)))
    expect_error(sampleControlElements(x, i.target, i.control, oversample = 100))
    expect_warning(s4 <- sampleControlElements(list(x, runif(length(x))),
                                               i.target, i.control, nbins = c(100,100)),
                   regexp = "Too few control elements")

    expect_equal(length(s1), length(i.target))
    expect_equal(s1, s2)
    expect_equal(s1, s3)
    expect_true(all(s1 %in% i.control))
    expect_equal(quantile(x[s1]), quantile(x[i.target]), tolerance = 0.1)
    expect_error(sampleControlElements(x, idxTarget = -1))
})
