context("plotGeneRegion")

test_that("prepareGTF fails with incorrect input", {
    gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf",
                           package = "swissknife")
    
    expect_error(prepareGTF())
    expect_error(prepareGTF(""))
    expect_error(prepareGTF("file.gtf"))
    expect_error(prepareGTF(gtffile, transcriptIdColumn = "nonexisting"))
    expect_error(prepareGTF(gtffile, geneIdColumn = "nonexisting"))
})

test_that("prepareGTF works", {
    gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf",
                           package = "swissknife")
    gtf <- prepareGTF(gtffile, transcriptIdColumn = "transcript_id",
                      geneIdColumn = "gene_id", geneSymbolColumn = "gene_name")
    expect_true("transcript" %in% colnames(S4Vectors::mcols(gtf)))
    expect_true("gene" %in% colnames(S4Vectors::mcols(gtf)))
    expect_true("gene_name" %in% colnames(S4Vectors::mcols(gtf)))
    expect_true(all(gtf$type == "exon"))
    
    expect_message(prepareGTF(gtffile, transcriptIdColumn = "transcript_id",
                              geneIdColumn = "gene_id", 
                              geneSymbolColumn = "missing"), 
                   "Using the 'gene' column as gene symbols")
})

test_that("plotGeneRegion fails with incorrect input", {
    gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf",
                           package = "swissknife")
    bwfile <- system.file("extdata/plotGeneRegion/mnase_mm10.bw",
                          package = "swissknife")
    names(bwfile) <- "bw1"
    
    expect_error(plotGeneRegion())
    expect_error(plotGeneRegion(gtf = 1))
    expect_error(plotGeneRegion(gtf = c("file1", "file2")))
    expect_error(plotGeneRegion(granges = 1))
    expect_error(plotGeneRegion(granges = gtffile))
    expect_error(plotGeneRegion(gtf = gtffile, chr = 1))
    expect_error(plotGeneRegion(gtf = gtffile, chr = c("chr1", "chr2")))
    expect_error(plotGeneRegion(gtf = gtffile, chr = "chr10", start = "1"))
    expect_error(plotGeneRegion(gtf = gtffile, chr = "chr10", start = 1, end = "2"))
    expect_error(plotGeneRegion(gtf = gtffile, chr = "chr10", start = c(1, 2)))
    expect_error(plotGeneRegion(gtf = gtffile, chr = "chr10", start = 1, end = c(1, 2)))
    expect_error(plotGeneRegion(gtf = gtffile, chr = "chr10", start = 200, end = 100))
    expect_error(plotGeneRegion(gtf = gtffile, showgene = 1))
    expect_error(plotGeneRegion(gtf = gtffile, showgene = c("gene1", "gene2")))
    expect_error(plotGeneRegion(gtf = gtffile, bigwigFiles = 1))
    expect_error(plotGeneRegion(gtf = gtffile, bigwigFiles = "file.bw"))
    expect_error(plotGeneRegion(gtf = gtffile, bigwigCond = 1))
    expect_error(plotGeneRegion(gtf = gtffile, bigwigCond = "c1"))
    expect_error(plotGeneRegion(gtf = gtffile, bigwigFiles = bwfile, bigwigCond = "c1"))
    expect_error(plotGeneRegion(gtf = gtffile, bigwigFiles = bwfile, bigwigCond = structure("c1", names = "n1")))
    expect_error(plotGeneRegion(gtf = gtffile, geneTrackTitle = 1))
    expect_error(plotGeneRegion(gtf = gtffile, geneTrackTitle = c("t1", "t2")))
    expect_error(plotGeneRegion(bigwigFiles = bwfile))
    expect_error(plotGeneRegion(gtf = gtffile))
    expect_error(plotGeneRegion(gtf = gtffile, lowerPadding = "1"))
    expect_error(plotGeneRegion(gtf = gtffile, upperPadding = "1"))
    expect_error(plotGeneRegion(gtf = gtffile, lowerPadding = c(1, 2)))
    expect_error(plotGeneRegion(gtf = gtffile, upperPadding = c(1, 2)))
})

test_that("plotGeneRegion works", {
    gtffile <- system.file("extdata/plotGeneRegion/mm10_ensembl98.gtf",
                           package = "swissknife")
    bwfile <- system.file("extdata/plotGeneRegion/mnase_mm10.bw",
                          package = "swissknife")
    names(bwfile) <- "bw1"
    bwcond <- c("c1")
    names(bwcond) <- "bw1"
    
    ## GTF file, no bigwig
    a1 <- plotGeneRegion(gtf = gtffile, showgene = "Bclaf1")
    expect_is(a1, "list")
    
    ## GTF file, no bigwig, region
    a11 <- plotGeneRegion(gtf = gtffile, showgene = "", chr = "chr10",
                          start = 20310000, end = 20340000)
    expect_is(a11, "list")
    
    ## GRanges, no bigwig
    b <- prepareGTF(gtffile)
    a2 <- plotGeneRegion(gtf = "", granges = b, showgene = "Bclaf1")
    expect_is(a2, "list")
    expect_identical(a1, a2)
    
    ## Both GTF and GRanges, no bigwig
    a3 <- plotGeneRegion(gtf = gtffile, granges = b, showgene = "Bclaf1")
    expect_is(a3, "list")
    expect_identical(a1, a3)
    
    ## No annotation, bigwig without grouping
    a4 <- plotGeneRegion(bigwigFiles = bwfile, chr = "chr10", start = 20310000, end = 20340000)
    expect_is(a4, "list")
    
    ## No annotation, bigwig with grouping
    a5 <- plotGeneRegion(bigwigFiles = bwfile, bigwigCond = bwcond, 
                         chr = "chr10", start = 20310000, end = 20340000)
    expect_is(a5, "list")
    
    ## GTF file, bigwig, gene provided
    a6 <- plotGeneRegion(gtf = gtffile, showgene = "Bclaf1", 
                         bigwigFiles = bwfile, bigwigCond = bwcond)
    expect_is(a6, "list")
    
    ## GTF file, bigwig, gene and coordinates provided
    a7 <- plotGeneRegion(gtf = gtffile, showgene = "Bclaf1", 
                         bigwigFiles = bwfile, bigwigCond = bwcond, 
                         chr = "chr10", start = 20310000, end = 20340000)
    expect_is(a7, "list")
    expect_equal(a6, a7)
})