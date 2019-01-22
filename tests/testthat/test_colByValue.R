context("colByValue")

test_that("colByValue() works properly", {

    # arguments
    expect_error(colByValue(x = "error"))
    expect_error(colByValue(x = 1:2, rng = "error"))
    expect_error(colByValue(x = 1:2, rng = 1:3))
    expect_error(colByValue(x = 1:2, col = "error"))
    expect_error(colByValue(x = 1:2, NA.col = "error"))
    expect_error(colByValue(x = 1:2, alpha = 256))
    expect_warning(colByValue(x = 1:2, col = "#10101020", alpha = 128))

    # return value
    res1 <- colByValue(x = 1:10)
    res2 <- colByValue(x = 1:10, rng = 4:5)
    res3 <- colByValue(x = 1:2, col = c("white","black"), alpha = 128)
    expect_is(res1, "character")
    expect_length(res1, 10L)
    expect_length(res2, 10L)
    expect_length(unique(res2), 2L)
    expect_identical(res3, c("#FFFFFF80", "#00000080"))
})
