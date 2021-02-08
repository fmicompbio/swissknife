test_that("getGenomicTiles() works properly", {
    ## simple test data
    genomefile <- system.file("extdata", "getMappableRegions", "hg19sub.fa", package = "swissknife")
    genomechrlen <- Biostrings::fasta.seqlengths(genomefile)
    genomename <- "BSgenome.HSapiens.QuasR.hg19sub"
    
    tss <- GenomicRanges::GRanges("chr1", IRanges(start = c(1000, 1505, 2001),
                                                  width = 1, names = paste0("t", 1:3)),
                                  strand = c("+", "-", "*"))
    cgi <- GenomicRanges::GRanges("chr2", IRanges(start = c(750, 3750), width = 370))
    tssL <- as(tss, "GRangesList")
    cgiL <- as(cgi, "GRangesList")

    expect_error(getGenomicTiles("error"))
    expect_error(getGenomicTiles(genomename))
    expect_error(getGenomicTiles(100000))
    expect_error(getGenomicTiles(genomechrlen, addSeqComp = FALSE))
    expect_warning(getGenomicTiles(genomechrlen, 1000, addSeqComp = TRUE))
    
    expect_s4_class(gr0 <- getGenomicTiles(c("chr1" = 1999), 500, addSeqComp = FALSE), "GRanges")
    expect_length(gr0, 3L)
    expect_s4_class(gr1 <- getGenomicTiles(genomechrlen, 500, addSeqComp = FALSE), "GRanges")
    expect_s4_class(gr2 <- getGenomicTiles(genomefile, 500, addSeqComp = FALSE), "GRanges")
    expect_true(all(GenomicRanges::width(gr1) == 500L))
    expect_length(gr1, 190)
    expect_identical(ncol(S4Vectors::mcols(gr1)), 0L)
    expect_identical(gr1, gr2)
    expect_s4_class(gr3 <- getGenomicTiles(genomefile, 500,
                                           hasOverlap = list(CGI1 = cgi),
                                           fracOverlap = list(CGI2 = cgi),
                                           numOverlap = list(TSS1 = tss),
                                           nearest = list(TSS2 = tss)),
                    "GRanges")
    expect_identical(ncol(S4Vectors::mcols(gr3)), 8L)
    df3 <- S4Vectors::mcols(gr3)
    expect_identical(colnames(df3),
                     c("percGC", "CpGoe", "CGI1.hasOverlap", "CGI2.fracOverlap",
                       "TSS1.numOverlapWithin", "TSS1.numOverlapAny",
                       "TSS2.nearestName", "TSS2.nearestDistance"))
    expect_true(all(df3$percGC >= 0 && df3$percGC <= 100))
    expect_true(all(df3$CpGoe > 0 && df3$CpGoe < 2))
    expect_identical(sum(df3$CGI1.hasOverlap), 4L)
    expect_true(all(df3$CGI1.hasOverlap == (df3$CGI2.fracOverlap > 0)))
    expect_identical(df3$CGI2.fracOverlap[df3$CGI1.hasOverlap],
                     rep(c(251 / 500, 119 / 500), 2))
    expect_identical(sum(df3$TSS1.numOverlapWithin), 3L)
    expect_identical(df3$TSS1.numOverlapWithin, df3$TSS1.numOverlapAny)
    expect_identical(which(as.logical(df3$TSS1.numOverlapWithin)), c(2L, 4L, 5L))
    expect_true(all(is.na(df3$TSS2.nearestName) == is.na(df3$TSS2.nearestDistance)))
    expect_identical(df3[1:5, "TSS2.nearestName"], rep(names(tss), c(3,1,1)))
    expect_identical(df3[1:5, "TSS2.nearestDistance"], c(499L, 0L, 0L, 0L, 0L))
    tss1 <- tss
    GenomicRanges::width(tss1) <- 100
    expect_s4_class(gr4 <- getGenomicTiles(genomechrlen, 500,
                                           numOverlap = list(TSS = tss1),
                                           addSeqComp = FALSE),
                    "GRanges")
    expect_true(all(gr4$TSS.numOverlapWithin <= gr4$TSS.numOverlapAny))
    expect_identical(gr4$TSS.numOverlapWithin[1:5], rep(c(0L, 1L), c(3, 2)))
    expect_identical(gr4$TSS.numOverlapAny[1:5], rep(c(0L, 1L), c(1, 4)))
    expect_s4_class(gr5 <- getGenomicTiles(genomefile, 500,
                                           hasOverlap = list(CGI1 = cgiL),
                                           fracOverlap = list(CGI2 = cgiL),
                                           numOverlap = list(TSS1 = tssL),
                                           nearest = list(TSS2 = tssL)),
                    "GRanges")
    expect_equal(gr3, gr5)
    
    ## extended test data (also test with a BSgenome)
    skip_if_not_installed("QuasR")
    skip_if_not_installed("BSgenome")

    library("BSgenome")

    # temporary R library
    rlibdir <- tempfile(pattern = "Rlib")
    dir.create(rlibdir)
    .libPaths(rlibdir)

    # install BSgenome.HSapiens.QuasR.hg19sub into temporary library
    bsgPkg <- system.file("extdata", "BSgenome.HSapiens.QuasR.hg19sub_0.1.0.tar.gz", package = "QuasR")
    suppressMessages({
        utils::install.packages(pkgs = bsgPkg, lib = rlibdir, repos = NULL,
                                type = "source", INSTALL_opts = "--no-test-load")
    })
    
    library(genomename, character.only = TRUE)
    genomeobj <- get(genomename)
    
    expect_error(getGenomicTiles(genomeobj))
    expect_s4_class(gr6 <- getGenomicTiles(genomename, 500,
                                           hasOverlap = list(CGI1 = cgi),
                                           fracOverlap = list(CGI2 = cgi),
                                           numOverlap = list(TSS1 = tss),
                                           nearest = list(TSS2 = tss)),
                    "GRanges")
    expect_s4_class(gr7 <- getGenomicTiles(genomeobj, 500,
                                           hasOverlap = list(CGI1 = cgi),
                                           fracOverlap = list(CGI2 = cgi),
                                           numOverlap = list(TSS1 = tss),
                                           nearest = list(TSS2 = tss)),
                    "GRanges")
    expect_equal(gr3, gr6)
    expect_equal(gr3, gr7)
})
