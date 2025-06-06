library(testthat)

## -------------------------------------------------------------------------- ##
## Checks, .assertScalar
## -------------------------------------------------------------------------- ##
test_that(".assertScalar works", {
    expect_error(.assertScalar(1, type = TRUE))
    expect_error(.assertScalar(1, type = 1))
    expect_error(.assertScalar(1, type = c("numeric", "character")))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = 1))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = c(TRUE, FALSE)))
    
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_true(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = 0))
    expect_true(.assertScalar(-1, type = "numeric", validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", validValues = c(-2, 0)))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_error(.assertScalar(NA, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 1)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 3:4))
    expect_true(.assertScalar(NULL, type = "numeric", allowNULL = TRUE))
    expect_error(.assertScalar(NULL, type = "numeric", allowNULL = FALSE))
    expect_error(.assertScalar(1, type = "character"))
    expect_error(.assertScalar("x", type = "numeric"))
    expect_error(.assertScalar(FALSE, type = "character"))
    expect_error(.assertScalar(c(1, 2), type = "numeric"))
    test <- "text"
    expect_error(.assertScalar(x = test, type = "numeric"),
                 ".test. must be of class .numeric.")
    expect_error(.assertScalar(x = list(a = 1)$a, type = "logical"),
                 "list(a = 1)$a", fixed = TRUE)
    tmp <- matrix(1:4, ncol = 2)
    expect_error(.assertScalar(x = tmp[, 1], type = "logical"),
                 "tmp[, 1]", fixed = TRUE)
})

## -------------------------------------------------------------------------- ##
## Checks, .assertVector
## -------------------------------------------------------------------------- ##
test_that(".assertVector works", {
    expect_error(.assertVector(1, type = TRUE))
    expect_error(.assertVector(1, type = 1))
    expect_error(.assertVector(1, type = c("numeric", "character")))
    expect_error(.assertVector(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertVector(1, type = "numeric", allowNULL = 1))
    expect_error(.assertVector(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertVector(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertVector(1, type = "numeric", allowNULL = c(TRUE, FALSE)))
    expect_error(.assertVector(1, type = "numeric", len = TRUE))
    expect_error(.assertVector(1, type = "numeric", len = "rng"))
    expect_error(.assertVector(1, type = "numeric", len = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngLen = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngLen = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1:3))
    
    expect_true(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 1.5)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_error(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 2)))
    expect_error(.assertVector(c(1, 2), type = "numeric", len = 1))
    expect_true(.assertVector(c(1, 2), type = "numeric", len = 2))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngLen = c(3, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(2, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(1, 2)))
    expect_error(.assertVector(c("a", "b"), type = "character", validValues = c("A", "B")))
    expect_true(.assertVector(LETTERS[1:2], type = "character", validValues = LETTERS))
    test <- "text"
    expect_error(.assertVector(x = test, type = "numeric"),
                 ".test. must be of class .numeric.")
})

## -------------------------------------------------------------------------- ##
## Checks, .assertPackagesAvailable
## -------------------------------------------------------------------------- ##
test_that(".assertPackagesAvailable works", {
    testfunc <- function(...) .assertPackagesAvailable(...)
    expect_error(testfunc(1L))
    expect_error(testfunc("test", "error"))
    expect_error(testfunc("test", c(TRUE, FALSE)))

    expect_true(testfunc("base"))
    expect_true(testfunc("githubuser/base"))
    expect_true(testfunc(c("base", "methods")))
    expect_error(testfunc(c("error", "error2")), "BiocManager")
    expect_error(testfunc("error1", suggestInstallation = FALSE),
                 "installed[.]")
    rm(testfunc)
})

## -------------------------------------------------------------------------- ##
## Checks, .message
## -------------------------------------------------------------------------- ##
test_that(".message works", {
    verbose <- TRUE
    expect_type(suppressMessages(.message("message")), "character")
    expect_type(suppressMessages(.message("message", noTimer = TRUE)), "character")
    verbose <- FALSE
    expect_null(suppressMessages(.message("message")))
})
