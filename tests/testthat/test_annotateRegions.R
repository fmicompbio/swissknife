test_that("annotateRegions() works properly", {
    ## simple test data
    x <- GenomicRanges::GRanges(rep(c("chr1", "chr2", "chr3"), each = 2),
                                IRanges(start = c(300, 1700, 800, 3200, 1700, 3000),
                                        width = 1000),
                                strand = c("+", "-", "*", "*", "+", "-"),
                                score = 1:6)
    tss <- GenomicRanges::GRanges("chr1", IRanges(start = c(1000, 1505, 2001),
                                                  width = 1, names = paste0("t", 1:3)),
                                  strand = c("-", "+", "*"))
    cgi <- GenomicRanges::GRanges("chr2", IRanges(start = c(750, 3750), width = 370))
    tssL <- as(tss, "GRangesList")
    cgiL <- as(cgi, "GRangesList")

    expect_error(annotateRegions(x = "error"))
    expect_error(annotateRegions(x, hasOverlap = "error"))
    expect_error(annotateRegions(x, fracOverlap = "error"))
    expect_error(annotateRegions(x, numOverlap = "error"))
    expect_error(annotateRegions(x, nearest = "error"))
    expect_error(annotateRegions(x, ignore.strand = "error"))

    expect_s4_class(x1 <- annotateRegions(x), "GRanges")
    expect_identical(x, x1)
    expect_s4_class(x2 <- annotateRegions(x,
                                          hasOverlap = list(CGI1 = cgi),
                                          fracOverlap = list(CGI2 = cgi),
                                          numOverlap = list(TSS1 = tss, CGI3 = cgi),
                                          nearest = list(TSS2 = tss),
                                          ignore.strand = TRUE),
                    "GRanges")
    expect_s4_class(x3 <- annotateRegions(x,
                                          hasOverlap = list(CGI1 = cgi),
                                          fracOverlap = list(CGI2 = cgi),
                                          numOverlap = list(TSS1 = tss, CGI3 = cgi),
                                          nearest = list(TSS2 = tss),
                                          ignore.strand = FALSE),
                    "GRanges")
    expect_s4_class(x4 <- annotateRegions(x,
                                          hasOverlap = list(CGI1 = cgiL),
                                          fracOverlap = list(CGI2 = cgiL),
                                          numOverlap = list(TSS1 = tssL, CGI3 = cgiL),
                                          nearest = list(TSS2 = tssL),
                                          ignore.strand = TRUE),
                    "GRanges")
    expect_s4_class(x5 <- annotateRegions(x,
                                          hasOverlap = list(CGI1 = cgiL),
                                          fracOverlap = list(CGI2 = cgiL),
                                          numOverlap = list(TSS1 = tssL, CGI3 = cgiL),
                                          nearest = list(TSS2 = tssL),
                                          ignore.strand = FALSE),
                    "GRanges")
    expect_identical(x2, x4)
    expect_identical(x3, x5)
    expect_false(identical(x2, x3))
    expect_length(x2, length(x))
    expect_length(x3, length(x))
    expect_identical(ranges(x2), ranges(x))
    expect_identical(ranges(x3), ranges(x))
    expect_identical(x2$score, x$score)
    expect_identical(x3$score, x$score)
    expect_identical(mcols(x2),
                     S4Vectors::DataFrame(
                         score = 1:6, 
                         CGI1.hasOverlap = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
                         CGI2.fracOverlap = c(0, 0, 0.32, 0.37, 0, 0),
                         TSS1.numOverlapWithin = c(1L, 1L, 0L, 0L, 0L, 0L),
                         TSS1.numOverlapAny = c(1L, 1L, 0L, 0L, 0L, 0L),
                         CGI3.numOverlapWithin = c(0L, 0L, 0L, 1L, 0L, 0L),
                         CGI3.numOverlapAny = c(0L, 0L, 1L, 1L, 0L, 0L),
                         TSS2.nearestName = c("t1", "t3", NA, NA, NA, NA),
                         TSS2.nearestDistance = c(0L, 0L, NA, NA, NA, NA)))
    expect_identical(mcols(x3),
                     S4Vectors::DataFrame(
                         score = 1:6, 
                         CGI1.hasOverlap = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
                         CGI2.fracOverlap = c(0, 0, 0.32, 0.37, 0, 0),
                         TSS1.numOverlapWithin = c(0L, 1L, 0L, 0L, 0L, 0L),
                         TSS1.numOverlapAny = c(0L, 1L, 0L, 0L, 0L, 0L),
                         CGI3.numOverlapWithin = c(0L, 0L, 0L, 1L, 0L, 0L),
                         CGI3.numOverlapAny = c(0L, 0L, 1L, 1L, 0L, 0L),
                         TSS2.nearestName = c("t2", "t3", NA, NA, NA, NA),
                         TSS2.nearestDistance = c(205L, 0L, NA, NA, NA, NA)))
})
