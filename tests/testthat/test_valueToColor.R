context("valueToColor")

test_that("valueToColor() works properly", {

    # arguments
    expect_error(valueToColor(x = "error"))
    expect_error(valueToColor(x = 1:2, rng = "error"))
    expect_error(valueToColor(x = 1:2, rng = 1:3))
    expect_error(valueToColor(x = 1:2, col = "error"))
    expect_error(valueToColor(x = 1:2, NA.col = "error"))
    expect_error(valueToColor(x = 1:2, alpha = 256))
    expect_warning(valueToColor(x = 1:2, col = "#10101020", alpha = 128))

    # return value
    res1 <- valueToColor(x = 1:10)
    res2 <- valueToColor(x = 1:10, rng = 4:5)
    res3 <- valueToColor(x = 1:2, col = c("white","black"), alpha = 128)
    expect_is(res1, "character")
    expect_length(res1, 10L)
    expect_length(res2, 10L)
    expect_length(unique(res2), 2L)
    expect_identical(res3, c("#FFFFFF80", "#00000080"))
})
