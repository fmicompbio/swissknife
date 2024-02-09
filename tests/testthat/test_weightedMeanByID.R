library(SummarizedExperiment)

set.seed(123)
meansRows <- sample(1:100, 10, replace = TRUE)
dat <- unlist(lapply(meansRows, function(m) {
    rnorm(n = 5, mean = m, sd = 0.1*m)
}))
ma <- matrix(dat, nrow = 10, ncol = 5, byrow = TRUE)
IDs <- data.frame(ID = sample(c("A", "B", "C", "D"), size = 10, replace = TRUE))
Groups <- data.frame(group = c("Y","Y", "Z", "Z", "Z"))

mockSE <- SummarizedExperiment(assays = list(promoterenr = ma),
                     rowData = IDs,
                     colData = Groups)


test_that("weightedMeanByID only takes appropriate input", {
    expect_error(weightedMeanByID(SummarizedExperiment::assay(mockSE), "promoterenr", idCol = "ID", log2Transformed = FALSE), "Provide data to aggregate in a SummarizedExperiment")
    expect_error(weightedMeanByID(mockSE, "counts", idCol = "ID", log2Transformed = FALSE), "The assay name must be an existing assay of the SummarizedExperiment object")
    expect_error(weightedMeanByID(mockSE, ma, idCol = "ID", log2Transformed = FALSE), "Provide the assay name as a character")
    expect_error(weightedMeanByID(mockSE, "promoterenr", log2Transformed = FALSE), "idCol must be an existing column name in the rowData of the SummarizedExperiment object")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = IDs, log2Transformed = FALSE), "Provide the name of the rowData column as a character")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = "ID", log2Transformed = "log10"), "log2Transformed must be a logical indicating if data is log2 transformed")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = "ID", groupCol = Groups, log2Transformed = FALSE), "Provide the name of the colData column as a character")
    expect_error(weightedMeanByID(mockSE, "promoterenr", idCol = "ID", groupCol = "fantastic", log2Transformed = FALSE), "groupCol must be an existing column name in the colData of the SummarizedExperiment object")
    
    alt_SE <- mockSE
    assay(alt_SE)[c(4,8)] <- NA
    expect_warning(weightedMeanByID(alt_SE, "promoterenr", idCol = "ID", log2Transformed = FALSE), "The assay contains NAs, remaining values will be aggregated.")
})


test_that("weightedMeanByID works properly", {
    res <- weightedMeanByID(mockSE, "promoterenr", idCol = "ID", log2Transformed = FALSE)
    expect_identical(length(unique(IDs$ID)), length(rownames(res)))
    expect_true(all(unique(IDs$ID) %in% rownames(res)))
    expect_true(all(unique(Groups$group) %in% colnames(res)))
    
    lapply(list("Y", "Z"), function(gc) {
        lapply(unique(IDs$ID), function(l) {
            sel <- which(rowData(mockSE)$ID == l)
            means <- lapply(sel, function(s) {
                sum(assay(mockSE)[s, colData(mockSE)$group == gc])/sum(colData(mockSE)$group == gc)
            })
            squares <- lapply(means, function(m) {
                m^2
            })
            expect_equal(res[l,gc], sum(unlist(squares))/sum(unlist(means)))
        })
    })
    
    expect_equal(ncol(res), length(unique(colData(mockSE)$group)))
    
    loggedSE <- mockSE
    assay(loggedSE) <- log2(assay(loggedSE))
    expect_equal(weightedMeanByID(loggedSE, "promoterenr", idCol = "ID", log2Transformed = TRUE), log2(res))
    expect_is(res, "data.frame")
    
    # test that result is NA if all input values are NA
    na_SE <- mockSE
    assay(na_SE)[rowData(na_SE)$ID == "A",] <- NA
    res_na <- weightedMeanByID(na_SE, "promoterenr", idCol = "ID", log2Transformed = FALSE)
    expect_true(all(is.na(res_na["A",])))
    
    # test that result is 0 if all input values are 0
    assay(na_SE)[rowData(na_SE)$ID == "A",] <- 0
    res_null <- weightedMeanByID(na_SE, "promoterenr", idCol = "ID", log2Transformed = FALSE)
    expect_true(all(res_null["A",] == 0))
    
    # test that this still works if we have only one group
    oneSE <- mockSE[, colData(mockSE)$group == "Y"]
    res_one <- weightedMeanByID(oneSE, "promoterenr", idCol = "ID", log2Transformed = FALSE)
    expect_identical(res$Y, res_one$Y)
})
