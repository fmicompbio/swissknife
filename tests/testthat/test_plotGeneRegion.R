context("plotGeneRegion")

test_that("plotGeneRegion fails with incorrect input", {
    gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf",
                           package = "swissknife")
    bwfile <- system.file("extdata/plotGeneRegion/mnase_mm10.bw",
                          package = "swissknife")
    names(bwfile) <- "bw1"
    
    expect_error(plotGeneRegion())
    expect_error(plotGeneRegion(gtf = 1))
    expect_error(plotGeneRegion(granges = 1))
    expect_error(plotGeneRegion(granges = gtffile))
})