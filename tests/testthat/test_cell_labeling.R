## create data
## ... synthetic SingleCellExperiment
set.seed(42)
mu <- ceiling(runif(n = 2000, min = 0, max = 100))
counts <- do.call(rbind, lapply(mu, function(x){rpois(300, lambda = x)}))
counts <- counts + 1
i <- sample(x = 1:nrow(counts), size = 100)
j <- sample(x = 1:ncol(counts), size = 100)
counts[i, j] <- counts[i, j] + sample(5:10, length(i), replace = TRUE)
colnames(counts) <- paste0("c", seq.int(ncol(counts)))
rownames(counts) <- paste0("g", seq.int(nrow(counts)))
## ... with three cell types
celltypes <- rep(c("tA","tB","tC"), each = 100)
celltypegenes <- split(sample(rownames(counts), size = 300),
                       celltypes)
markergenes <- lapply(celltypegenes, "[", 1:10)
for (ct in unique(celltypes)) {
    counts[celltypegenes[[ct]], celltypes == ct] <- 2 * counts[celltypegenes[[ct]], celltypes == ct]
}
sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts))
## ... and log2 normalized counts
libsizes <- colSums(counts)
SingleCellExperiment::sizeFactors(sce) <- libsizes/mean(libsizes)

SingleCellExperiment::logcounts(sce) <- log2(sweep(as.matrix(SingleCellExperiment::counts(sce)),
                                                   2, SingleCellExperiment::sizeFactors(sce), "/") + 1)


test_that("normGenesetExpression() works properly", {
    ## these should fail
    expect_error(normGenesetExpression("not_a_sce", NULL))
    expect_error(normGenesetExpression(sce, c("1","2","3")))
    expect_error(normGenesetExpression(sce, c("not","known","genes")))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), expr_values = FALSE))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = "error"))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = 1:3))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = -1))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), nbins = "error"))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), nbins = 1:3))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), nbins = -1))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = 20, subset.row = c(TRUE, FALSE)))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = 20, subset.row = c(-1, 10)))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = 20, subset.row = c("g99999")))
    expect_error(normGenesetExpression(sce, c("g1","g2","g3"), R = 20, BPPARAM = "error"))
    
    ## run function
    set.seed(7)
    selgns <- sample(rownames(sce), size = 20)
    res1 <- normGenesetExpression(sce = sce, genes = selgns,
                                  expr_values = "logcounts", R = 40,
                                  nbins = 100, subset.row = rownames(sce),
                                  BPPARAM = BiocParallel::SerialParam())
    res2 <- normGenesetExpression(sce = sce, genes = selgns,
                                  expr_values = 2, R = 40,
                                  nbins = 100, subset.row = seq.int(nrow(sce)),
                                  BPPARAM = BiocParallel::MulticoreParam(workers = 4L))
    res3 <- normGenesetExpression(sce = sce[unique(c(paste0("g",1:100),selgns))], genes = selgns,
                                  expr_values = 2, R = 2, nbins = 100, 
                                  subset.row = rep(TRUE, 120L))
    
    ## tests
    expect_type(res1, "double")
    expect_length(res1, ncol(sce))
    expect_gt(cor(res1, res2), 0.9)
    expect_identical(names(res1), names(res2))
    expect_identical(names(res1), names(res3))
    expect_identical(names(res1), colnames(sce))
})

test_that("labelCells() works properly", {
    ## these should fail
    expect_error(labelCells("not_a_sce", NULL))
    expect_error(labelCells(sce, c("1","2","3")))
    expect_error(labelCells(sce, list("1","2","3")))
    expect_error(labelCells(sce, list(tA = 1:2)))
    expect_error(labelCells(sce, list(tA = 1:2, tB = 3)))
    expect_error(labelCells(sce, list(tA = 1:2, tB = 3:4)))
    expect_error(labelCells(sce, list(tA = c("1","2"), tB = c("3","4"))))
    expect_error(labelCells(sce, markergenes, expr_values = FALSE))
    expect_error(labelCells(sce, markergenes, fraction_topscoring = "low"))
    expect_error(labelCells(sce, markergenes, fraction_topscoring = c(.01, .05)))
    expect_error(labelCells(sce, markergenes, fraction_topscoring = 2.5))
    expect_error(labelCells(sce, markergenes, normGenesetExpressionParams = list(TRUE)))
    expect_error(labelCells(sce, markergenes, aggregateReferenceParams = list(TRUE)))
    expect_error(labelCells(sce, markergenes, SingleRParams = list(TRUE)))

    ## run function
    res1 <- labelCells(sce = sce, markergenes = markergenes,
                       expr_values = "logcounts",
                       normGenesetExpressionParams = list(R = 40),
                       BPPARAM = BiocParallel::SerialParam())
    expect_warning({
    res2 <- labelCells(sce = sce, markergenes = markergenes,
                       expr_values = 2,
                       normGenesetExpressionParams = list(sce = sce, R = 40),
                       aggregateReferenceParams = list(power = 0.5, BPPARAM = NULL),
                       SingleRParams = list(BPPARAM = NULL),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 4L))
    })
    
    ## tests
    expect_type(res1, "list")
    expect_length(res1, 3L)
    expect_named(res1, c("cells", "refs", "labels"))
    expect_true(all(diag(table(res1$labels$pruned.labels, celltypes)) > 95))
    expect_true(all(diag(table(res2$labels$pruned.labels, celltypes)) > 95))
})
