#' @title Aggregate different rows assigned to the same ID by calculating a 
#' weighted mean 
#'
#' @description First row means are calculated to summarize across replicates 
#' identified by the groupCol in the colData. Then different row means that are 
#' assigned to the same feature ID given by the idCol in the rowData are 
#' summarized by calculating a weighted mean. This weighted mean is the sum of 
#' the squared row means divided by the sum of the row means. If all row means 
#' are 0, they remain 0 in the output.
#' 
#' @author Fiona Ross
#' 
#' @param SE a \code{SummarizedExperiment} object that contains an assay with
#' values to be aggregated, a colData column that assigns samples to their group
#' and a rowData column with IDs to indicate which rows to combine. 
#' @param assay the name of the assay in the SummarizedExperiment object that 
#' should be aggregated.
#' @param idCol the column name in the rowData of the SummarizedExperiment 
#' indicating the feature ID.
#' @param groupCol the column name in the colData of the SummarizedExperiment 
#' indicating which columns belong to the same group and should be averaged as
#' replicates, before the weighted mean is calculated across rows.
#' @param log2Transformed a \code{logical} indicating whether values in the 
#' assay are log2 transformed. If log2Transformed is TRUE, an exponential
#' transformation will be applied before aggregating the values and another log 
#' transformation afterwards. 
#' 
#' @return The output is a \code{data.frame} with one column for each of the 
#' unique names in the groupCol and one row for each of the unique IDs in the
#' idCol. The row and column names are the respective unique values. The entries
#' represent the weighted means for each unique feature ID. If all the input 
#' values were NA, the aggregated value is also NA, while for all zero, the 
#' output remains zero. If log2Transformed is true the output will be log2
#' transformed again.
#' 
#' @examples
#' set.seed(123)
#' meansRows <- sample(1:100, 10, replace = TRUE)
#' dat <- unlist(lapply(meansRows, function(m) {
#'     rnorm(n = 5, mean = m, sd = 0.1*m)
#' }))
#' ma <- matrix(dat, nrow = 10, ncol = 5, byrow = TRUE)
#' IDs <- data.frame(ID = sample(c("A", "B", "C", "D"), size = 10, replace = TRUE))
#' Groups <- data.frame(group = c("Y","Y", "Z", "Z", "Z"))
#' mockSE <- SummarizedExperiment(assays = list(promoterenr = ma),
#'                               rowData = IDs,
#'                               colData = Groups)
#' weightedMeanByID(mockSE, "counts", idCol = "ID", log2Transformed = FALSE)                                
#' 
#' @importFrom stats aggregate weighted.mean
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom dplyr mutate across where
#' @export
weightedMeanByID <- function(SE, assay, idCol = "GENEID", groupCol = "group", log2Transformed = TRUE) {
    .assertPackagesAvailable("stats")
    if (!is(SE, "SummarizedExperiment")) {
        stop("Provide data to aggregate in a SummarizedExperiment")
    }
    if (!is.character(assay)) {
        stop("Provide the assay name as a character")
    }
    if (!is.character(idCol)) {
        stop("Provide the name of the rowData column as a character")
    }
    if (!is.character(groupCol)) {
        stop("Provide the name of the colData column as a character")
    }
    if (!is.logical(log2Transformed)) {
        stop("log2Transformed must be a logical indicating if data is log2 transformed")
    }
    if (!idCol %in% colnames(SummarizedExperiment::rowData(SE))) {
        stop("idCol must be an existing column name in the rowData of the SummarizedExperiment object")
    }
    if (!groupCol %in% colnames(SummarizedExperiment::colData(SE))) {
        stop("groupCol must be an existing column name in the colData of the SummarizedExperiment object")
    }
    if (!assay %in% names(SummarizedExperiment::assays(SE))) {
        stop("The assay name must be an existing assay of the SummarizedExperiment object")
    }
    if (any(is.na(SummarizedExperiment::assay(SE)))) {
        warning("The assay contains NAs, remaining values will be aggregated.")
    }
    
    if (log2Transformed) {
        ma <- 2^SummarizedExperiment::assay(SE, assay)
    } else {
        ma <- SummarizedExperiment::assay(SE, assay)
    }
    
    weightedMeanByID <- do.call(cbind, lapply(unique(SummarizedExperiment::colData(SE)[, groupCol]), function(g) {
        sel <- colData(SE)[, groupCol] == g
        ma <- ma[, sel, drop = FALSE]
        promoterenr <- as.data.frame(rowMeans(ma, na.rm = TRUE))
        promoterenr$gene_id <- SummarizedExperiment::rowData(SE)[, idCol]
        
        weightedMeanByID <- stats::aggregate(
            . ~ gene_id,
            data = promoterenr,
            FUN = function(x) {
                if (sum(x) != 0) {
                    stats::weighted.mean(x, x)
                } else if (sum(x) == 0) {
                    0
                }
            })
        
        if (log2Transformed) {
            weightedMeanByID <- weightedMeanByID %>%
                dplyr::mutate(across(where(is.numeric), log2))
        }
        
        rownames(weightedMeanByID) <- weightedMeanByID$gene_id
        weightedMeanByID$gene_id <- NULL
        colnames(weightedMeanByID) <- c("weightedMean")
        
        if (any(!unique(promoterenr$gene_id) %in% rownames(weightedMeanByID))) {
            append <- unique(promoterenr$gene_id)[!unique(promoterenr$gene_id) %in% rownames(weightedMeanByID)]
            append <- data.frame(weightedMean = rep(NA, length(append)), row.names = append)
            weightedMeanByID <- rbind(weightedMeanByID, append)
        }
        
        return(weightedMeanByID)
    }))
    colnames(weightedMeanByID) <- unique(SummarizedExperiment::colData(SE)[, groupCol])
    return(weightedMeanByID)
}
