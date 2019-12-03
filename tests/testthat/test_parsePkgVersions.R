context("parsePkgVersions")

test_that("parsePkgVersions() works properly", {
    infiles <- list.files(system.file("extdata", "parsePkgVersions",
                                      package = "swissknife"),
                          full.names = TRUE)
    badfiles <- tempfile(fileext = c(".bad", ".Rout"))
    writeLines(letters, badfiles[2])

    ## pre-flight checks
    expect_error(parsePkgVersions(1L))
    expect_error(parsePkgVersions(badfiles))
    expect_warning(res <- parsePkgVersions(badfiles[2]))
    expect_identical(res, list())

    ## Trying to read a non-existent file should give a warning and an empty data.frame
    expect_is(res <- parsePkgVersions(infiles), "list")
    expect_length(res, 2L)
    expect_match(names(res), "^R version [0-9]\\.[0-9]\\.[0-9]")
    expect_identical(unname(lengths(res)), c(2L, 2L))
    expect_identical(lengths(res[[1]]), c(files = 1L, packages = 23L))
    expect_identical(lengths(res[[2]]), c(files = 1L, packages = 17L))
    expect_identical(order(names(res)), 1:2)
    expect_identical(order(res[[1]]$packages), seq_along(res[[1]]$packages))
    expect_true(all(sapply(res, "[[", "files") %in% infiles))
    expect_identical(res[[2]]$packages[16], "tidyr_1.0.0")
})
