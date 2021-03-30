test_that(".getDistMat() works properly", {

    ## create data
    x <- 1:1000
    x2 <- x
    names(x2) <- paste0("x", 1:1000)
    y <- sort(rnorm(n = 1000, mean = 0, sd = 5)) + runif(n = 1000, min = 0, max = 10)
    lo <- loess(y ~ x)
    x_curve <- seq(range(x)[1], range(x)[2], length.out=100)
    
    ## run function
    euclMat <- .getDistMat(loessModel=lo, x=x, y=y, x_curve=x_curve)
    euclMat2 <- .getDistMat(loessModel=lo, x=x, y=y, x_curve=x_curve, method="manhattan")
    euclMat3 <- .getDistMat(loessModel=lo, x=x2, y=y, x_curve=x_curve)
    
    ## tests
    expect_true(length(x)==nrow(euclMat))
    expect_true(length(x_curve)==ncol(euclMat))
    expect_error(.getDistMat(x=x, y=y, x_curve=x_curve))
    expect_error(.getDistMat(loessModel=lo, y=y, x_curve=x_curve))
    expect_error(.getDistMat(loessModel=lo, x=x, x_curve=x_curve))
    expect_error(.getDistMat(loessModel=lo, x=x, y=y))
    expect_true(!is.null(euclMat2))
    expect_true(all(rownames(euclMat3)==names(x2)))
    
})


test_that("selVarGenes() works properly", {
    
    ## create data
    mu <- ceiling(runif(n = 2000, min = 0, max = 100))
    counts <- do.call(rbind, lapply(mu, function(x){rpois(1000, lambda = x)}))
    counts <- counts + 1
    i <- sample(x = 1:nrow(counts), size = 500)
    j <- sample(x = 1:ncol(counts), size = 500)
    counts[i, j] <- counts[i, j] + sample(5:10, length(i), replace = TRUE)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts))
    libsizes <- colSums(counts)
    SingleCellExperiment::sizeFactors(sce) <- libsizes/mean(libsizes)
    
    SingleCellExperiment::logcounts(sce) <- log2(sweep(as.matrix(SingleCellExperiment::counts(sce)), 2, SingleCellExperiment::sizeFactors(sce), "/") + 1)
    
    
    ## run function
    sel <- selVarGenes(data = sce, assay.type = "counts")
    sel2 <- selVarGenes(data = 2^SingleCellExperiment::logcounts(sce)-1)
    sel3 <- selVarGenes(data = sce, assay.type = "logcounts")
    
    ## tests
    expect_error(selVarGenes())
    expect_error(selVarGenes(data = mu))
    expect_error(selVarGenes(data = sce, assay.type = "normCounts"))
    expect_error(selVarGenes(data = SingleCellExperiment::SingleCellExperiment(assays = list(counts=counts))))
    expect_error(selVarGenes(data = SingleCellExperiment::SingleCellExperiment(assays = list())))
    expect_error(selVarGenes(data = SingleCellExperiment::SingleCellExperiment(assays = list(counts=counts)), assay.type="logcounts"))
    expect_true(all(sel$varGenes%in%rownames(sel$geneInfo)))
    expect_true(!is.null(sel$varGenes))
    expect_true(all(sel$varGenes==sel2$varGenes))
    expect_true(all(sel$varGenes==sel3$varGenes))
})


test_that("plotSelVarGenes() works properly", {
    
    ## create data
    mu <- ceiling(runif(n = 2000, min = 0, max = 100))
    counts <- do.call(rbind, lapply(mu, function(x){rpois(1000, lambda = x)}))
    counts <- counts + 1
    i <- sample(x = 1:nrow(counts), size = 500)
    j <- sample(x = 1:ncol(counts), size = 500)
    counts[i, j] <- counts[i, j] + sample(5:10, length(i), replace = TRUE)
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=counts))
    libsizes <- colSums(counts)
    SingleCellExperiment::sizeFactors(sce) <- libsizes/mean(libsizes)
    
    ## run selVarGenes
    sel <- selVarGenes(data = sce)
    
    ## tests
    expect_true(plotSelVarGenes(sel))
    expect_true(plotSelVarGenes(sel, colByBin=TRUE))
    expect_error(plotSelVarGenes())
    expect_error(plotSelVarGenes(mu))
    expect_error(plotSelVarGenes(list(mu=mu)))
    
})



