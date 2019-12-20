context("selVarGenes")

test_that("getDistMat() works properly", {
     
     ## create data
     x <- 1:1000
     x2 <- x
     names(x2) <- paste0("x", 1:1000)
     y <- sort(rnorm(n = 1000, mean = 0, sd = 5)) + runif(n = 1000, min = 0, max = 10)
     lo <- loess(y ~ x)
     x_curve <- seq(range(x)[1], range(x)[2], length.out=100)
     
     ## run function
     euclMat <- getDistMat(loessModel=lo, x=x, y=y, x_curve=x_curve)
     euclMat2 <- getDistMat(loessModel=lo, x=x, y=y, x_curve=x_curve, method="manhattan")
     euclMat3 <- getDistMat(loessModel=lo, x=x2, y=y, x_curve=x_curve)
     
     ## tests
     expect_true(length(x)==nrow(euclMat))
     expect_true(length(x_curve)==ncol(euclMat))
     expect_error(getDistMat(x=x, y=y, x_curve=x_curve))
     expect_error(getDistMat(loessModel=lo, y=y, x_curve=x_curve))
     expect_error(getDistMat(loessModel=lo, x=x, x_curve=x_curve))
     expect_error(getDistMat(loessModel=lo, x=x, y=y))
     expect_true(!is.null(euclMat2))
     expect_true(all(rownames(euclMat3)==names(x2)))

})

