test_that("calcPhasogram(), estimateNRL() and calcAndCountDist() work properly", {
    ## create phasograms using MNase data from chr10:20000000-21000000
    bamf <- system.file("extdata", "phasograms", "mnase_mm10.bam", package = "swissknife")
    pg1 <- calcPhasogram(bamf)
    pg2 <- calcPhasogram(bamf, rmdup = FALSE)
    pg3 <- calcPhasogram(c(bamf,bamf))
    pg4 <- calcPhasogram(bamf, regions = GenomicRanges::GRanges("chr10", IRanges::IRanges(start = 19000000, end = 23000000)))
    pg5 <- calcPhasogram(bamf, regions = GenomicRanges::GRanges("chr11", IRanges::IRanges(start = 19000000, end = 23000000)))

    ## check invalid arguments
    expect_warning(estimateNRL(rep(0,3000L))) # all zeros
    expect_warning(estimateNRL(pg1, usePeaks = 1:100)) # too many usePeaks
    expect_error(calcAndCountDist(1:3, 3:1, numeric(3))) # not sorted ascendingly

    ## check if the plotting function runs
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)
    expect_true(is.list(plotPhasogram(x = pg1)))
    expect_true(is.list(plotPhasogram(x = pg1, verbosePlot = TRUE)))
    dev.off()
    unlink(tf)

    ## check expected results
    expect_true(is.numeric(pg1))
    expect_length(pg1, 3000L)
    expect_true(all(pg2 > pg1))
    expect_equal(pg1, pg3)
    expect_equal(pg1, pg4)
    expect_true(sum(pg5) == 0L)

    nrl <- estimateNRL(pg1, usePeaks = 1:4)
    expect_true(is.list(nrl))
    expect_equal(nrl$nrl, 186.0)
    expect_length(nrl$nrl.CI95, 2L)

    expect_equal(calcAndCountDist(c(1,2,4),c(1,3,5),numeric(4)), c(2,1,1,1))
    expect_equal(calcAndCountDist(c(1,3,5),c(1,2,4),numeric(4)), c(2,0,1,0))
})
