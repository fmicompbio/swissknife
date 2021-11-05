test_that("getInsertSizeDistFromBam() works properly", {
    # bam file
    bamf <- system.file("extdata", "getInsertSizeDistFromBam", "atac_mm10.bam",
                        package = "swissknife")
    
    bh <- Rsamtools::scanBamHeader(bamf)[[1]]$targets
    bh <- bh[setdiff(names(bh), c("chrM", "chrX", "chrY"))]
    chrsReg <- GenomicRanges::GRanges(names(bh),
                                      IRanges::IRanges(start = 1L,
                                                       width = bh,
                                                       names = names(bh)))
    
    
    # argument checks
    expect_error(getInsertSizeDistFromBam(fname = "error"), "failed to open")
    expect_error(getInsertSizeDistFromBam(fname = bamf, regions = "error"), "regions")
    expect_error(getInsertSizeDistFromBam(fname = bamf, nmax = "error"), "nmax")
    expect_error(getInsertSizeDistFromBam(fname = bamf, isizemax = 1:2), "isizemax")
    expect_error(getInsertSizeDistFromBam(fname = bamf, exclude = NA), "exclude")
    
    # correct results
    expect_is(res1 <- getInsertSizeDistFromBam(bamf), "integer")
    expect_is(res2 <- getInsertSizeDistFromBam(c(bamf, bamf)), "integer")
    expect_is(expect_message(res3 <- getInsertSizeDistFromBam(bamf, regions = chrsReg)), "integer")
    expect_is(res4 <- getInsertSizeDistFromBam(bamf, nmax = 12361), "integer")
    expect_is(expect_message(res5 <- getInsertSizeDistFromBam(bamf, isizemax = 400)), "integer")

    expect_identical(sum(res1), 13577L)

    expect_length(res1, 800L)
    expect_length(res2, 800L)
    expect_length(res3, 800L)
    expect_length(res4, 800L)
    expect_length(res5, 400L)

    expect_identical(2L * res1, res2)
    expect_identical(res1, res3)
    expect_identical(sum(res4), 12361L)
    expect_identical(res1[1:399], res5[1:399])
})