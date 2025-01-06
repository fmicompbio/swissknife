# prepare test data
genomefile <- system.file("extdata", "getMappableRegions", "hg19sub.fa", package = "swissknife")
chrs <- Biostrings::readDNAStringSet(genomefile)
indexdir <- tempfile()
indexpre <- "index"
indexname <- file.path(indexdir, indexpre)
suppressMessages({
    idx <- Rbowtie::bowtie_build(genomefile, indexdir)
})

test_that(".getChrlenFromBowtieIndex works properly", {
    expect_type(len <- .getChrlenFromBowtieIndex(indexname), "integer")
    expect_identical(len, structure(width(chrs), names = names(chrs)))
})

test_that(".writeWindowsToTempFile works properly", {
    chrs1 <- as.character(XVector::subseq(chrs[[1]], start = 1, width = 100))
    tf1 <- tempfile()
    expect_type(res1 <- .writeWindowsToTempFile(chr = chrs1, w = 10, fname = tf1), "character")
    expect_identical(tf1, res1)
    expect_identical(length(lns1 <- readLines(tf1)), 182L)
    expect_identical(lns1[c(1,2,181,182)],
                     c(">1", substr(chrs1, 1, 10),
                       ">91", substr(chrs1, 91, 100)))
    expect_type(res2 <- .writeWindowsToTempFile(chr = chrs1, w = 100), "character")
    expect_identical(length(lns2 <- readLines(res2)), 2L)
    expect_identical(lns2, c(">1", chrs1))
    unlink(c(res1, res2))
})

test_that(".writeWindowsToTempFileCPP works properly", {

    chrs1 <- as.character(XVector::subseq(chrs[[1]], start = 1, width = 100))
    expect_error(.writeWindowsToTempFileCPP(chr = chrs1, w = 10, fname = file.path(tempdir(), "error", "error.fa")))
    tf1 <- tempfile()
    expect_type(res1 <- .writeWindowsToTempFileCPP(chr = chrs1, w = 10, fname = tf1), "character")
    expect_identical(tf1, res1)
    expect_identical(length(lns1 <- readLines(tf1)), 182L)
    expect_identical(lns1[c(1,2,181,182)],
                     c(">1", substr(chrs1, 1, 10),
                       ">91", substr(chrs1, 91, 100)))

    expect_type(res2 <- .writeWindowsToTempFileCPP(chr = chrs1, w = 100, tempfile()), "character")
    expect_identical(length(lns2 <- readLines(res2)), 2L)
    expect_identical(lns2, c(">1", chrs1))

    res3 <- .writeWindowsToTempFile(chr = chrs1, w = 10, fname = tempfile())
    lns3 <- readLines(res3)
    expect_identical(lns1, lns3)

    unlink(c(res1, res2, res3))
})

test_that(".alignWindowsToGenome works properly", {
    skip_if_not_installed("Rbowtie")

    chrs1 <- as.character(XVector::subseq(chrs[[1]], start = 2600, width = 1000))
    tf1 <- .writeWindowsToTempFile(chr = chrs1, w = 50)

    tfs <- c(tempfile(), tempfile(), tempfile())
    expect_type(res1 <- .alignWindowsToGenome(tf1, indexname, m = 1, p = 2),
              "character")
    expect_named(res1, c("fmax", "fun", "fout"))
    expect_identical(length(readLines(res1["fmax"])), 590L)
    expect_type(res2 <- .alignWindowsToGenome(tf1, indexname, m = 10, p = 2,
                                            fmax = tfs[1], fun = tfs[2], fout = tfs[3]),
              "character")
    expect_identical(tfs, unname(res2))
    expect_identical(length(readLines(res2["fmax"])), 158L)
})

test_that("getMappableRegions() works properly", {
    readfile <- tempfile(fileext = ".fa")
    writeLines(c(">test", "ACGTACGTCATGCTGACTGACTGACGA"), readfile)

    expect_error(getMappableRegions("error", indexname, 50))
    expect_error(getMappableRegions(genomefile, file.path(indexname, indexpre), 50))
    expect_error(getMappableRegions(genomefile, file.path(indexdir, "error"), 50))
    expect_error(getMappableRegions(genomefile, FALSE, 50))
    expect_error(getMappableRegions(readfile, indexname, 50))
    expect_message(gr1 <- getMappableRegions(genomefile, indexname, 50, quiet = FALSE))
    expect_s4_class(gr1, "GRanges")
    expect_length(gr1, 27L)
    expect_identical(sum(width(gr1)), 94090L)

    skip_if_not_installed("QuasR")
    skip_if_not_installed("BSgenome")
    skip_if_not_installed("withr")

    # copy sample data
    file.copy(system.file("extdata", package = "QuasR"), tempdir(), recursive = TRUE)

    # -------------- temporary R library
    withr::local_temp_libpaths()
    withr::local_envvar(c(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)))

    # create BSgenome and Rbowtie index
    clObj <- parallel::makeCluster(2L)
    bsgPkg <- file.path(tempdir(), "extdata", "BSgenome.HSapiens.QuasR.hg19sub_0.1.0.tar.gz")
    suppressMessages({
        utils::install.packages(pkgs = bsgPkg, repos = NULL, quiet = TRUE,
                                type = "source", INSTALL_opts = "--no-test-load")
    })
    samplefile <- file.path(tempdir(), "extdata", "samples_chip_single.txt")
    genomePkg <- "BSgenome.HSapiens.QuasR.hg19sub"
    proj <- QuasR::qAlign(samplefile, genomePkg, clObj = clObj, lib.loc = .libPaths()[1])
    parallel::stopCluster(cl = clObj)

    # test getMappableRegions
    expect_s4_class(gr2 <- getMappableRegions(genomePkg, indexname, 50), "GRanges")
    expect_identical(gr1, gr2)
    expect_s4_class(gr3 <- getMappableRegions(genomePkg, paste0(genomePkg, ".Rbowtie"), 50), "GRanges")
    expect_identical(gr1, gr3)
    expect_s4_class(gr4 <- getMappableRegions(get(genomePkg), indexname, 50), "GRanges")
    expect_identical(gr1, gr4)

    # -------------- restore original library paths
    withr::deferred_run()

    unlink(readfile)
})
