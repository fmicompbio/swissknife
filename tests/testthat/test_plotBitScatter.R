test_that("plotBitScatter() runs", {
    x <- rnorm(100)
    y <- rnorm(100)

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_null(plotBitScatter(x = x, y = y))
    expect_null(plotBitScatter(x = x, y = y, col = "gray"))

    dev.off()
    unlink(tf)
})
