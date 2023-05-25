test_that("loadExampleData() works properly", {

    # arguments
    expect_error(loadExampleData(FALSE), "name")
    expect_error(loadExampleData(c("mycars", "mycars_summary")), "name")
    expect_error(loadExampleData("error"), "must be one of")
    expect_error(loadExampleData("mycars", "error"), "envir")
    expect_error(loadExampleData("mycars", globalenv(), "error"), "verbose")

    # listing
    expect_output(df0 <- loadExampleData("list"))
    expect_is(df0, "data.frame")
    expect_identical(colnames(df0),
                     c("name", "created", "description", "is_collection"))

    # return value
    expect_message(df1 <- loadExampleData("mycars", verbose = TRUE))
    expect_is(df1, "data.frame")
    expect_true(exists("mycars", envir = globalenv()))
    expect_identical(dim(get("mycars", envir = globalenv())), c(32L, 13L))
    rm("mycars", envir = globalenv())

    expect_message(df2 <- loadExampleData("mycars_summary", verbose = TRUE))
    expect_is(df2, "data.frame")
    expect_true(exists("mycars_summary", envir = globalenv()))
    expect_identical(dim(get("mycars_summary", envir = globalenv())), c(3L, 3L))
    rm("mycars_summary", envir = globalenv())

    df3 <- loadExampleData("2023-05", verbose = FALSE)
    expect_is(df3, "data.frame")
    expect_identical(dim(df3), c(2L, 4L))
    expect_identical(rbind(df1, df2), df3)
    expect_true(exists("mycars", envir = globalenv()))
    expect_true(exists("mycars_summary", envir = globalenv()))
    rm("mycars", "mycars_summary", envir = globalenv())

    expect_message(df4 <- loadExampleData("new", verbose = TRUE))
    expect_is(df4, "data.frame")
    expect_identical(df3, df4)
    expect_true(exists("mycars", envir = globalenv()))
    expect_true(exists("mycars_summary", envir = globalenv()))
    rm("mycars", "mycars_summary", envir = globalenv())

    # environment
    myenv <- new.env()
    df1 <- loadExampleData("mycars", envir = myenv, verbose = FALSE)
    expect_true(!exists("mycars", envir = globalenv()))
    expect_true(exists("mycars", envir = myenv))
    expect_identical(dim(get("mycars", envir = myenv)), c(32L, 13L))
    rm("mycars", envir = myenv)
    rm(myenv)
})