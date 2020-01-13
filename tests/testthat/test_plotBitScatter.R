context("plotBitScatter")

test_that("plotBitScatter() runs", {
    
    options(bitmapType = "cairo")
    
    x <- rnorm(100)
    y <- rnorm(100)

    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)

    expect_null(plotBitScatter(x = x, y = y))
    expect_null(plotBitScatter(x = x, y = y, col = "gray"))
    expect_null(plotBitScatter(cbind(x, y)))

    dev.off()
    unlink(tf)
})
