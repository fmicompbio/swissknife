# calculate measures of tissue-specific gene expression


#' @title Calculate gene-expression specificity scores.
#' 
#' @description Calculate expression specificity scores for genes that quantify
#'    specific expression of a gene in groups of samples (e.g. from different tissues).
#'
#' @author Michael Stadler
#'
#' @param x Expression data, either a \code{matrix} with expression values for
#'   genes (rows) in each sample (columns), or a \code{SummarizedExperiment} or
#'   \code{SingleCellExperiment} object containing such expression data in one
#'   of the assays (selected by \code{expr_values}).
#' @param method \code{character} scalar selecting the type of expression
#'   specificity score to be calculated. One of: "tau", "TSI", "counts".
#'   See "Details" for method-specific information.
#' @param group \code{character} or \code{factor} of length \code{ncol(x)} that
#'   groups the measurements into clusters or tissues, for which expression
#'   specificity scores are to be calculated. If \code{NULL} (the default), each
#'   column of \code{x} is assumed to be its own group. If multiple columns belong
#'   to the same group, these columns are first aggregated (averaged) before
#'   score calculations.
#' @param thresh \code{numeric} scalar defining the expression threshold. Values
#'   greater than this threshold are interpreted as expressed (used only for
#'   some of the methods, see "Details").
#' @param expr_values Integer scalar or string indicating which assay of
#'   \code{x} contains the expression values, for \code{x} of type
#'   \code{SummarizedExperiment} or \code{SingleCellExperiment}. Ignored if
#'   \code{x} is a \code{matrix}.
#' @param na.rm Logical scalar. If \code{TRUE}, \code{NA} values are excluded in
#'   the calculations.
#' 
#' @return A \code{numeric} vector of length \code{nrow(x)} with gene scores.
#' 
#' @references For a review of tissue-specificity scores, see:
#' "A benchmark of gene expression tissue-specificity metrics"
#' Nadezda Kryuchkova-Mostacci and Marc Robinson-Rechavi
#' Brief Bioinform. 2017 Mar; 18(2): 205–214.
#' doi: 10.1093/bib/bbw008, PMCID: PMC5444245, PMID: 26891983
#' 
#' @examples 
#' x <- rbind(g1 = runif(5),
#'            g2 = c(1, 0, 0, 0, 0),
#'            g3 = c(.6, .1, .1, .1, .1))
#' specificityScore(x)
#' specificityScore(x, method = "TSI")
#' specificityScore(x, method = "counts", thresh = 0.5)
#' 
#' @importFrom matrixStats rowAvgsPerColSet rowSums2 rowMaxs
#' @importFrom methods .valueClassTest
#' 
#' @export
#' @docType methods
#' @rdname specificityScore-methods
setGeneric(name = "specificityScore",
           def = function(x, method = c("tau", "TSI", "counts"), group = NULL,
                          thresh = 0, expr_values = "logcounts", na.rm = FALSE) {
               standardGeneric("specificityScore")
           },
           valueClass = "numeric")


#' @rdname specificityScore-methods
#' @aliases specificityScore,matrix-method
setMethod("specificityScore",
          c("matrix"),
          function(x, method = c("tau", "TSI", "counts"), group = NULL,
                   thresh = 0, expr_values = "logcounts", na.rm = FALSE) {
              # pre-flight checks
              method <- match.arg(method)
              stopifnot(exprs = {
                  is.matrix(x)
                  is.numeric(thresh)
                  length(thresh) == 1L
                  is.logical(na.rm)
                  length(na.rm) == 1L
              })
              
              # aggregate observations by group
              if (!is.null(group)) {
                  stopifnot(exprs = {
                      is.character(group) || is.factor(group) || is.numeric(group)
                      length(group) == ncol(x)
                  })
                  message("aggregating observations by group")
                  if (is.character(group) || is.numeric(group)) {
                      group <- factor(as.character(group),
                                      levels = unique(as.character(group)))
                  }
                  ml <- max(table(group))
                  S <- do.call(cbind, lapply(split(seq.int(ncol(x)), group),
                                             function(y) {
                                                 c(y, rep(NA, ml - length(y)))
                                             }))
                  x <- rowAvgsPerColSet(X = x, S = S)
              }
              
              # calculate expression specificity score
              message("calculating '", method, "' for ", nrow(x), " genes and ",
                      ncol(x), " groups")
              scores <- switch(method,
                  tau = rowSums2(1 - (x / rowMaxs(x, na.rm = na.rm)), na.rm = na.rm) / (ncol(x) - 1),
                  TSI = rowMaxs(x, na.rm = na.rm) / rowSums2(x, na.rm = na.rm),
                  counts = rowSums2(x > thresh, na.rm = na.rm)
              )
              return(scores)
          })


#' @rdname specificityScore-methods
#' @aliases specificityScore,SummarizedExperiment-method
setMethod("specificityScore",
          c("SummarizedExperiment"),
          function(x, method = c("tau", "TSI", "counts"), group = NULL,
                   thresh = 0, expr_values = "logcounts", na.rm = FALSE) {
              .assertPackagesAvailable("SummarizedExperiment")
              ## check expr_values
              if (is.numeric(expr_values)) {
                  stopifnot(exprs = {
                      length(expr_values) == 1L
                      expr_values <= length(SummarizedExperiment::assays(x))
                  })
              } else if (is.character(expr_values)) {
                  stopifnot(exprs = {
                      length(expr_values) == 1L
                      expr_values %in% SummarizedExperiment::assayNames(x)
                  })
              } else {
                  stop("'expr_values' is not a valid value for use in assay()")
              }
              ## extract expression matrix
              xx <- as.matrix(SummarizedExperiment::assay(x, expr_values))
              ## call matrix method
              specificityScore(xx, method, group, thresh, expr_values, na.rm)
          })
