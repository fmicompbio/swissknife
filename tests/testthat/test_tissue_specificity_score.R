context("tissue specificity score")

## create data
## ... synthetic matrix
set.seed(42)
x <- rbind(matrix(runif(n = 90, min = 0, max = 100), ncol = 5),
           rep(c(0,100), c(4, 1)),
           rep(c(100,0), c(1, 4)))
x.NA <- x; x.NA[sample(length(x), size = 5)] <- NA
x2 <- x[, rep(seq.int(ncol(x)), each = 2)]
se <- SummarizedExperiment::SummarizedExperiment(list(e = x))
sce <- as(se, "SingleCellExperiment")


test_that("specificityScore() works properly", {
    ## these should fail
    expect_error(specificityScore("not_a_valid_x"))
    expect_error(specificityScore(x, method = "not_a_valid_method"))
    expect_error(specificityScore(x, "counts", group = "error"))
    expect_error(specificityScore(x, "counts", thresh = "error"))
    expect_error(specificityScore(x, "counts", thresh = 1:2))
    expect_error(specificityScore(x, "counts", na.rm = c(TRUE,FALSE)))
    expect_error(specificityScore(x, "counts", na.rm = "error"))
    expect_error(specificityScore(se, expr_values = 3))
    expect_error(specificityScore(se, expr_values = "non_existing"))
    expect_error(specificityScore(se, expr_values = TRUE))
    expect_error(specificityScore(sce, expr_values = 3))
    expect_error(specificityScore(sce, expr_values = "non_existing"))
    expect_error(specificityScore(sce, expr_values = TRUE))
    
    ## run function
    ## ... tau
    res.tau1 <- specificityScore(x, method = "tau")
    res.tau2 <- specificityScore(se, method = "tau", expr_values = 1)
    res.tau3 <- specificityScore(se, method = "tau", expr_values = "e")
    res.tau4 <- specificityScore(sce, method = "tau", expr_values = 1)
    res.tau5 <- specificityScore(sce, method = "tau", expr_values = "e")
    res.tau6 <- specificityScore(x.NA, method = "tau")
    res.tau7 <- specificityScore(x.NA, method = "tau", na.rm = TRUE)
    res.tau8 <- specificityScore(x2, method = "tau",
                                 group = rep(seq.int(ncol(x)), each = 2))
    ## ... TSI
    res.TSI1 <- specificityScore(x, method = "TSI")
    res.TSI2 <- specificityScore(x2, method = "TSI",
                                 group = rep(seq.int(ncol(x)), each = 2))
    ## ... counts
    res.counts1 <- specificityScore(x, method = "counts", thresh = min(x))
    res.counts2 <- specificityScore(x, method = "counts", thresh = 50)
    res.counts3 <- specificityScore(x2, method = "counts", thresh = 50,
                                    group = rep(seq.int(ncol(x)), each = 2))
    
    ## tests
    ## ... tau
    expect_is(res.tau1, "numeric")
    expect_length(res.tau1, nrow(x))
    expect_true(all(res.tau1 >= 0 & res.tau1 <= 1))
    expect_identical(res.tau1, res.tau2)
    expect_identical(res.tau1, res.tau3)
    expect_identical(res.tau1, res.tau4)
    expect_identical(res.tau1, res.tau5)
    expect_true(any(is.na(res.tau6)))
    expect_identical(res.tau1[!is.na(res.tau6)], res.tau6[!is.na(res.tau6)])
    expect_true(all(!is.na(res.tau7)))
    expect_identical(res.tau1[!is.na(res.tau6)], res.tau7[!is.na(res.tau6)])
    expect_true(all(res.tau1[is.na(res.tau6)] != res.tau7[is.na(res.tau6)]))
    expect_identical(res.tau1, res.tau8)
    expect_equal(res.tau1, rowSums(1 - (x / apply(x, 1, max))) / (ncol(x) - 1))
    expect_identical(res.tau1[19:20], c(1, 1))
    ## ... TSI
    expect_is(res.TSI1, "numeric")
    expect_length(res.TSI1, nrow(x))
    expect_true(all(res.TSI1 >= 0 & res.TSI1 <= 1))
    expect_identical(res.TSI1, res.TSI2)
    expect_equal(res.TSI1, apply(x, 1, max) / rowSums(x))
    expect_identical(res.TSI1[19:20], c(1, 1))
    ## ... counts
    expect_is(res.counts1, "numeric")
    expect_length(res.counts1, nrow(x))
    expect_true(all(res.counts1 >= 0 & res.counts1 <= ncol(x)))
    expect_true(all(res.counts1 >= res.counts2))
    expect_identical(res.counts2, res.counts3)
    expect_equal(res.counts1, rowSums(x > 0))
    expect_equal(res.counts2, rowSums(x > 50))
    expect_identical(res.counts1, rep(c(ncol(x), 1), c(18, 2)))
})
