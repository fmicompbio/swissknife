test_that("sampleControlElements() works properly", {
    ## create synthetic data
    n.control <- 10000
    n.target <- 500
    i.control <- seq.int(n.control)
    i.target <- seq.int(n.target) + n.control
    set.seed(0)
    x <- c(runif(n.control, min = 0, max = 10),
           rnorm(n.target, mean = 5, sd = 1))

    ## check invalid arguments
    expect_error(sampleControlElements(letters, 1:10)) # non-list/non-numeric x
    expect_error(sampleControlElements(list(1:10, letters[1:10]), 1:5)) # list x with non-numeric elements
    expect_warning(sampleControlElements(list(x, runif(length(x))), # too many bins
                                         i.target, i.control, nbins = c(100,100)),
                   regexp = "Too few control elements")
    expect_error(sampleControlElements(x, i.target, i.control, oversample = 100))
    expect_error(sampleControlElements(x, idxTarget = -1))

    ## check expected results (univariate)
    set.seed(1)
    s1 <- sampleControlElements(x, idxTarget = i.target, idxControl = i.control)
    set.seed(1)
    s2 <- sampleControlElements(x, idxTarget = rep(c(FALSE,TRUE),c(n.control,n.target)), idxControl = NULL)
    set.seed(1)
    s3 <- sampleControlElements(list(x), idxTarget = rep(c(FALSE,TRUE),c(n.control,n.target)), idxControl = i.control)
    set.seed(1)
    s4 <- sampleControlElements(list(x), idxTarget = i.target, idxControl = rep(c(TRUE,FALSE),c(n.control,n.target)))
    expect_type(s1, "integer")
    expect_equal(length(s1), n.target)
    expect_equal(s1, s2)
    expect_equal(s1, s3)
    expect_equal(s1, s4)
    expect_true(all(s1 %in% i.control))
    expect_equal(quantile(x[s1]), quantile(x[i.target]), tolerance = 0.1)

    ## check expected results (multivariate)
    set.seed(2)
    s5 <- sampleControlElements(list(x, x), idxTarget = i.target, idxControl = i.control, nbins = 20)
    expect_type(s5, "integer")
    expect_equal(length(s5), n.target)
    expect_true(all(s5 %in% i.control))
    expect_equal(quantile(x[s5]), quantile(x[i.target]), tolerance = 0.1)
})
