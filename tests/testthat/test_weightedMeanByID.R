library(SummarizedExperiment)

set.seed(123)
meansRows <- sample(1:100, 10, replace = TRUE)
dat <- unlist(lapply(meansRows, function(m) {
    rnorm(n = 3, mean = m, sd = 0.1*m)
}))
ma <- matrix(dat, nrow = 10, ncol = 3, byrow = TRUE)
IDs <- data.frame(ID = sample(c("A", "B", "C", "D"), size = 10, replace = TRUE))

mockSE <- SummarizedExperiment(assays = list(promoterenr = ma),
                     rowData = IDs)


test_that("weightedMeanByID only takes appropriate input", {
    expect_error(weightedMeanByID(SummarizedExperiment::assay(mockSE), "promoterenr", idCol = "ID", log2Transformed = FALSE), "Provide data to aggregate in a Summarized Experiment")
    expect_error(weightedMeanByID(mockSE, "counts", idCol = "ID", log2Transformed = FALSE), "The assay name must be an existing assay of the Summarized Experiment object")
    expect_error(weightedMeanByID(mockSE, ma, idCol = "ID", log2Transformed = FALSE), "Provide the name of the Summarized Experiment assay to use for weighted mean calculation as a character")
    expect_error(weightedMeanByID(mockSE, "promoterenr", log2Transformed = FALSE), "The ID column must be an existing columnname in the rowData of the Summarized Experiment object")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = IDs, log2Transformed = FALSE), "Provide the name of the ID column over which the weighted mean is calculated as a character")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = "ID", log2Transformed = "log10"), "log2Transformed can only be TRUE or FALSE to indicate whether the data is log2 transformed")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = "ID", log2Transformed = TRUE, pseudocount = TRUE), "Provide the pseudocount as a numeric")
    
    alt_SE <- mockSE
    assay(alt_SE)[c(4,8)] <- NA
    expect_warning(weightedMeanByID(alt_SE, "promoterenr", idCol = "ID", log2Transformed = TRUE), "There are NAs in the assay to aggregate. rowMeans will be calculated over remaining values and used for the weighted mean.")
})


test_that("weightedMeanByID works properly", {
    res <- weightedMeanByID(mockSE, "promoterenr", idCol = "ID", log2Transformed = FALSE)
    expect_identical(length(unique(IDs$ID)), length(rownames(res)))
    
    lapply(unique(IDs$ID), function(l) {
        sel <- which(rowData(mockSE)$ID == l)
        means <- lapply(sel, function(s) {
            sum(assay(mockSE)[s,])/3 + 0.01
        })
        squares <- lapply(means, function(m) {
            m^2
        })
        expect_identical(res[l,], sum(unlist(squares))/sum(unlist(means)))
    })
    expect_equal(ncol(res), 1)
    loggedSE <- mockSE
    assay(loggedSE) <- log2(assay(loggedSE))
    expect_equal(weightedMeanByID(loggedSE, "promoterenr", idCol = "ID", log2Transformed = TRUE), log2(res))
    expect_is(res, "data.frame")
})
