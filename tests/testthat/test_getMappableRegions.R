context("getMappableRegions")

# prepare test data
library(Rbowtie)
library(Biostrings)
library(BSgenome)
genomefile <- system.file("extdata", "getMappableRegions", "hg19sub.fa", package = "swissknife")
chrs <- readDNAStringSet(genomefile)
indexdir <- tempfile()
indexpre <- "index"
indexname <- file.path(indexdir, indexpre)
idx <- bowtie_build(genomefile, indexdir)

test_that(".getChrlenFromBowtieIndex works properly", {
    expect_type(len <- .getChrlenFromBowtieIndex(indexname), "integer")
    expect_identical(len, structure(width(chrs), names = names(chrs)))
})

test_that(".writeWindowsToTempFile works properly", {
    chrs1 <- as.character(subseq(chrs[[1]], start = 1, width = 100))
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

test_that(".alignWindowsToGenome works properly", {
    chrs1 <- as.character(subseq(chrs[[1]], start = 2600, width = 1000))
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

    if (require(QuasR)) {

        # copy sample data
        file.copy(system.file("extdata", package = "QuasR"), tempdir(), recursive = TRUE)
        
        # temporary R library
        rlibdir <- tempfile(pattern = "Rlib")
        dir.create(rlibdir)
        .libPaths(rlibdir)
        # ... add rlibdir to R_LIBS for "R CMD INSTALL" and cluster nodes to find it
        oldRlibs <- Sys.getenv("R_LIBS")
        Sys.setenv(R_LIBS = paste(tools::file_path_as_absolute(rlibdir), oldRlibs,
                                  sep = .Platform$path.sep))
        
        # create cluster object
        clObj <- parallel::makeCluster(2L)
        
        # load QuasR on cluster nodes
        parallel::clusterEvalQ(cl = clObj, expr = library(QuasR))
        
        # install BSgenome.HSapiens.QuasR.hg19sub into temporary library
        bsgPkg <- file.path(tempdir(), "extdata", "BSgenome.HSapiens.QuasR.hg19sub_0.1.0.tar.gz")
        utils::install.packages(pkgs = bsgPkg, lib = rlibdir, repos = NULL,
                                type = "source", INSTALL_opts = "--no-test-load")
        
        # build Rbowtie index
        samplefile <- file.path(tempdir(), "extdata", "samples_chip_single.txt")
        genomePkg <- "BSgenome.HSapiens.QuasR.hg19sub"
        proj <- qAlign(samplefile, genomePkg, clObj = clObj, lib.loc = rlibdir)
            
        # test getMappableRegions
        expect_s4_class(gr2 <- getMappableRegions(genomePkg, indexname, 50),
                        "GRanges")
        expect_identical(gr1, gr2)
        expect_s4_class(gr3 <- getMappableRegions(genomePkg, paste0(genomePkg, ".Rbowtie"), 50),
                        "GRanges")
        expect_identical(gr1, gr3)
        expect_s4_class(gr4 <- getMappableRegions(get(genomePkg), indexname, 50),
                        "GRanges")
        expect_identical(gr1, gr4)
    }
    
    unlink(readfile)
})
