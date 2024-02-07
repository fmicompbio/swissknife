#' @title Aggregate different rows assigned to the same ID by calculating a 
#' weighted mean 
#'
#' @description After calculating row means to summarize across replicates, 
#' values of different rows that are assigned to the same ID are summarized by 
#' calculating a weighted mean. The weights are the values themselves.
#' 
#' @author Fiona Ross
#' 
#' @param SE \code{SummarizedExperiment} object that contains an assay with 
#' the values to be aggregated and a rowData column with IDs to indicate
#' which values of different rows to combine. 
#'   
#' @param assay \code{character} name of the assay in the SummarizedExperiment
#' object to use.
#' 
#' @param idCol \code{character} column name in the rowData of the
#' SummarizedExperiment indicating the ID assigned to each row.
#' 
#' @param log2Transformed \code{logical} indicating whether the data is log2
#' transformed. 
#' 
#' @param pseudocount \code{numeric} pseudocount to add to avoid dividing by 0.
#' 
#' @return \code{data.frame} with one column that is a numeric vector of the
#' weighted means per unique ID. If log2Transformed is true the output will be 
#' log2 transformed again. The rownames of the dataframe correspond to the 
#' unique IDs.
#' 
#' @importFrom stats aggregate weighted.mean
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom dplyr mutate across where
 
 
weightedMeanByID <- function(SE, assay, idCol = "GENEID", log2Transformed = TRUE, pseudocount = 0.01) {
    .assertPackagesAvailable("stats")
    if (!is(SE, "SummarizedExperiment")) {
        stop("Provide data to aggregate in a Summarized Experiment")
    }
    if (!is.character(assay)) {
        stop("Provide the name of the Summarized Experiment assay to use for weighted mean calculation as a character")
    }
    if (!is.character(idCol)) {
        stop("Provide the name of the ID column over which the weighted mean is calculated as a character")
    }
    if (!is.logical(log2Transformed)) {
        stop("log2Transformed can only be TRUE or FALSE to indicate whether the data is log2 transformed")
    }
    if (!idCol %in% colnames(SummarizedExperiment::rowData(SE))) {
        stop("The ID column must be an existing columnname in the rowData of the Summarized Experiment object")
    }
    if (!assay %in% names(SummarizedExperiment::assays(SE))) {
        stop("The assay name must be an existing assay of the Summarized Experiment object")
    }
    if (!is.numeric(pseudocount)) {
        stop("Provide the pseudocount as a numeric")
    }
    if (any(is.na(SummarizedExperiment::assay(SE)))) {
        warning("There are NAs in the assay to aggregate. rowMeans will be calculated over remaining values and used for the weighted mean.")
    }
    
    if (log2Transformed) {
        ma <- 2^SummarizedExperiment::assay(SE, assay)
    } else {
        ma <- SummarizedExperiment::assay(SE, assay)
    }
    
    promoterenr <- as.data.frame(rowMeans(ma, na.rm = TRUE))
    promoterenr$gene_id <- SummarizedExperiment::rowData(SE)[,idCol]
    
    weightedMeanByID <- stats::aggregate(
        . ~ gene_id,
        data = promoterenr,
        FUN = function(x)stats::weighted.mean(x + pseudocount, x + pseudocount))
    
    if (log2Transformed) {
        weightedMeanByID <- weightedMeanByID %>%
            dplyr::mutate(across(where(is.numeric), log2))
    }
    
    rownames(weightedMeanByID) <- weightedMeanByID$gene_id
    weightedMeanByID$gene_id <- NULL
    colnames(weightedMeanByID) <- c("weightedMean")
    
    return(weightedMeanByID)
}
