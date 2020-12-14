#' @title Calculate normalized expression of a gene set
#' 
#' @description Calculate normalized expression for a set of genes in each cell
#'   from a \code{SingleCellExperiment}, using random sets of similarly
#'   expressed genes as background to account for cell quality and
#'   sequencing depth.
#'
#' @author Michael Stadler
#'
#' @param sce \code{SingleCellExperiment} object.
#' @param genes \code{character} vector with the genes in the set. Must be a
#'   subset of \code{rownames(sce)}.
#' @param expr_values Integer scalar or string indicating which assay of
#'   \code{sce} contains the expression values.
#' @param R Integer scalar giving the number of random gene sets to sample
#'   for normalization.
#' @param subset.row Sample random genes only from these. If \code{NULL}
#'   (the default), the function will sample from all genes in \code{sce}.
#'   Alternatively, \code{subset.row} can be a logical, integer or character
#'   vector indicating the rows (genes) of \code{sce} to use for sampling.
#'   This allows for example to exclude highly variable genes from the sampling
#'   which are likely expressed only in certain cell types.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'   instance determining the parallel back-end to be used during evaluation.
#' 
#' @return A \code{numeric} vector with normalized gene set scores for each
#'   cell in \code{sce}.
#' 
#' @examples 
#' if (require(SingleCellExperiment)) {
#'     # get sce
#'     example(SingleCellExperiment, echo=FALSE)
#'     rownames(sce) <- paste0("g", seq.int(nrow(sce)))
#'     
#'     # calculate gene set expression scores
#'     markers <- c("g1", "g13", "g27")
#'     scores <- normGenesetExpression(sce, markers, R = 50)
#'     
#'     # compare expression of marker genes with scores
#'     plotdat <- cbind(scores, t(logcounts(sce)[markers, ]))
#'     cor(plotdat)
#'     pairs(plotdat)
#' }
#' 
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom matrixStats rowMeans2 rowSds
#' 
#' @export
normGenesetExpression <- function(sce,
                                  genes,
                                  expr_values = "logcounts",
                                  subset.row = NULL,
                                  R = 200,
                                  BPPARAM = SerialParam()) {
    # pre-flight checks
    stopifnot(exprs = {
        is(sce, "SingleCellExperiment")
        is.character(genes)
        length(genes) > 1L
        all(genes %in% rownames(sce))
        is.null(subset.row) || !any(is.na(subset.row))
        is.numeric(R)
        length(R) == 1L
        R > 0
        inherits(BPPARAM, "BiocParallelParam")
    })
    if (is.numeric(expr_values)) {
        stopifnot(exprs = {
            length(expr_values) == 1L
            expr_values <= length(assays(sce))
        })
    } else if (is.character(expr_values)) {
        stopifnot(exprs = {
            length(expr_values) == 1L
            expr_values %in% assayNames(sce)
        })
    } else {
        stop("'expr_values' is not a valid value for use in assay()")
    }

    # exclude genes if subset.row is given
    expr <- as.matrix(assay(sce, expr_values))
    i <- match(genes, rownames(sce))
    if (!is.null(subset.row)) {
        if (is.logical(subset.row) && length(subset.row) == nrow(sce)) {
            sel <- which(subset.row)
        } else if (is.numeric(subset.row) && min(subset.row) >= 0 && max(subset.row) <= nrow(sce)) {
            sel <- subset.row
        } else if (is.character(subset.row) && all(subset.row %in% rownames(sce))) {
            sel <- match(subset.row, rownames(sce))
        } else {
            stop("'subset.row' is not a valid selector of rows (genes) in 'sce'")
        }
        sel <- unique(c(i, sel)) # keep marker genes
        expr <- expr[sel, ]
        i <- match(genes, rownames(expr))
    }
    
    # calculate average expression of selected genes
    avgExpr <- rowMeans(expr)
    val.obs <- colSums(expr[i, , drop = FALSE])
    
    # calculate expression in random gene sets
    val.rand <- do.call(cbind, bplapply(seq_len(R), function(r) {
        irand <- sampleControlElements(x = avgExpr, idxTarget = i, nbins = 100)
        colSums(expr[irand, , drop = FALSE])
    }, BPPARAM = BPPARAM))
    
    # normalize observed gene set expression
    # z values: (obs - exp) / sqrt(Var[exp])
    mus <- rowMeans2(val.rand)
    sds <- rowSds(val.rand)
    val.obs.norm <- (val.obs - mus) / sds
    names(val.obs.norm) <- colnames(sce)
    
    return(val.obs.norm)
}


#' @title Assign labels to cells using known marker genes
#' 
#' @description Given marker gene sets for cell types, identify cells with
#'   high expression of the marker genes (positive examples), then use these
#'   cells to create a reference transcriptome profile for each cell type and
#'   identify additional cells of each type using \code{SingleR}. These
#'   marker genes should specifically expressed a single cell type, e.g.
#'   CD3 which is expressed by all T cell subtypes would not be suitable
#'   for specific T cell subtypes.
#'
#' @author Michael Stadler
#'
#' @param sce \code{SingleCellExperiment} object.
#' @param markergenes Named \code{list} of \code{character} vectors with the
#'   marker genes for each cell types. The marker genes must be a subset of
#'   \code{rownames(sce)}.
#' @param fraction_topscoring \code{numeric} vector of length 1 or the same
#'   length as \code{markergenes} giving the fraction(s) of top scoring cells
#'   for each cell type to pick to create the reference transcriptome profile.
#' @param expr_values Integer scalar or string indicating which assay of
#'   \code{sce} contains the expression values.
#' @param normGenesetExpressionParams \code{list} with additional parameters
#'   for \code{\link{normGenesetExpression}}.
#' @param aggregateReferenceParams \code{list} with additional parameters
#'   for \code{\link[SingleR]{aggregateReference}}.
#' @param SingleRParams \code{list} with additional parameters for
#'   \code{\link[SingleR]{SingleR}}.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'   instance determining the parallel back-end to be used during evaluation.
#' 
#' @return A \code{list} of three elements named \code{cells}, \code{refs} and \code{labels}.
#'   \code{cells} contains a \code{list} with the numerical indices of the top
#'   scoring cells for each cell type.
#'   \code{refs} contains the pseudo-bulk transcriptome profiles used as a 
#'   reference for label assignment, as returned by \code{\link[SingleR]{aggregateReference}}.
#'   \code{labels} contains a \code{\link[S4Vectors]{DataFrame}} with the
#'   annotation statistics for each cell (one cell per row), generated by
#'   \code{\link[SingleR]{SingleR}}.
#'   
#' @seealso \code{\link{normGenesetExpression}} used to calculate scores for
#'    marker gene sets; \code{\link[SingleR]{aggregateReference}} used to
#'    create reference profiles; \code{\link[SingleR]{SingleR}} used to assign
#'    labels to cells.
#' 
#' @examples 
#' # create SingleCellExperiment with cell-type specific genes
#' library(SingleCellExperiment)
#' n_types <- 3
#' n_per_type <- 30
#' n_cells <- n_types * n_per_type
#' n_genes <- 500
#' fraction_specific <- 0.1
#' n_specific <- round(n_genes * fraction_specific)
#' 
#' set.seed(42)
#' mu <- ceiling(runif(n = n_genes, min = 0, max = 30))
#' u <- do.call(rbind, lapply(mu, function(x) rpois(n_cells, lambda = x)))
#' rownames(u) <- paste0("g", seq.int(nrow(u)))
#' celltype.labels <- rep(paste0("t", seq.int(n_types)), each = n_per_type)
#' celltype.genes <- split(sample(rownames(u), size = n_types * n_specific),
#'                         rep(paste0("t", seq.int(n_types)), each = n_specific))
#' for (i in seq_along(celltype.genes)) {
#'     j <- celltype.genes[[i]]
#'     k <- celltype.labels == paste0("t", i)
#'     u[j, k] <- 2 * u[j, k]
#' }
#' v <- log2(u + 1)
#' sce <- SingleCellExperiment(assays=list(counts=u, logcounts=v))
#'
#' # define marker genes (subset of true cell-type-specific genes)
#' marker.genes <- lapply(celltype.genes, "[", 1:5)
#' marker.genes
#' 
#' # predict cell types
#' res <- labelCells(sce, marker.genes,
#'                   fraction_topscoring = 0.1,
#'                   normGenesetExpressionParams = list(R = 50))
#' 
#' # high-scoring cells used as references for each celltype
#' res$cells
#' 
#' # ... from these, pseudo-bulks were created:
#' res$refs
#' 
#' # ... and used to predict labels for all cells
#' res$labels$pruned.labels
#' 
#' # compare predicted to true cell types
#' table(true = celltype.labels, predicted = res$labels$pruned.labels)
#'       
#' @importFrom SingleR SingleR aggregateReference
#' @importFrom BiocParallel SerialParam bplapply
#' 
#' @export
labelCells <- function(sce,
                       markergenes,
                       fraction_topscoring = 0.01,
                       expr_values = "logcounts",
                       normGenesetExpressionParams = list(R = 200),
                       aggregateReferenceParams = list(power = 0.5),
                       SingleRParams = list(),
                       BPPARAM = SerialParam()) {
    ## pre-flight checks
    stopifnot(exprs = {
        # sce
        is(sce, "SingleCellExperiment")
        # markergenes
        is(markergenes, "list")
        !is.null(names(markergenes))
        length(markergenes) > 1L
        all(lengths(markergenes) > 1L)
        all(vapply(markergenes, is.character, logical(1)))
        all(unlist(markergenes) %in% rownames(sce))
        # fraction_topscoring
        is.numeric(fraction_topscoring)
        length(fraction_topscoring) == 1L || length(fraction_topscoring) == length(markergenes)
        all(fraction_topscoring > 0) && all(fraction_topscoring < 1)
        # normGenesetExpressionParams
        is(normGenesetExpressionParams, "list")
        length(normGenesetExpressionParams) == 0 || !is.null(names(normGenesetExpressionParams))
        # aggregateReferenceParams
        is(aggregateReferenceParams, "list")
        length(aggregateReferenceParams) == 0 || !is.null(names(aggregateReferenceParams))
        # SingleRParams
        is(SingleRParams, "list")
        length(SingleRParams) == 0 || !is.null(names(SingleRParams))
    })
    if (is.numeric(expr_values)) {
        stopifnot(exprs = {
            length(expr_values) == 1L
            expr_values <= length(assays(sce))
        })
    } else if (is.character(expr_values)) {
        stopifnot(exprs = {
            length(expr_values) == 1L
            expr_values %in% assayNames(sce)
        })
    } else {
        stop("'expr_values' is not a valid value for use in assay()")
    }

    ## calculate scores for marker gene sets
    todrop <- intersect(names(normGenesetExpressionParams),
                        c("sce", "genes", "expr_values", "BPPARAM"))
    if (length(todrop) > 0) {
        warning("ignoring the following user-provided elements in ",
                "'normGenesetExpressionParams': ",
                paste(todrop, collapse = ", "))
        normGenesetExpressionParams[todrop] <- NULL
    }
    scoreL <- lapply(markergenes, function(gns)
        do.call(normGenesetExpression,
                c(list(sce = sce, genes = gns,
                       expr_values = expr_values,
                       BPPARAM = BPPARAM),
                  normGenesetExpressionParams))
    )
    
    ## identify top scoring cells of each type and create references
    if (length(fraction_topscoring) == 1L)
        fraction_topscoring <- rep(fraction_topscoring, length(markergenes))
    topL <- lapply(seq_along(scoreL), function(i) {
        which(scoreL[[i]] > quantile(scoreL[[i]],
                                     probs = 1 - fraction_topscoring[i]))
    })
    names(topL) <- names(scoreL)
    topcells <- unlist(topL, use.names = FALSE)
    topcelltypes <- factor(rep(names(topL), lengths(topL)),
                           levels = names(markergenes))

    todrop <- intersect(names(aggregateReferenceParams),
                        c("ref", "labels", "BPPARAM"))
    if (length(todrop) > 0) {
        warning("ignoring the following user-provided elements in ",
                "'aggregateReferenceParams': ",
                paste(todrop, collapse = ", "))
        aggregateReferenceParams[todrop] <- NULL
    }
    refs <- do.call(aggregateReference,
                    c(list(ref = sce[, topcells],
                           labels = topcelltypes,
                           BPPARAM = BPPARAM),
                      aggregateReferenceParams))

    # transfer labels using SingleR
    todrop <- intersect(names(SingleRParams),
                        c("test", "ref", "labels", "BPPARAM"))
    if (length(todrop) > 0) {
        warning("ignoring the following user-provided elements in ",
                "'SingleRParams': ",
                paste(todrop, collapse = ", "))
        SingleRParams[todrop] <- NULL
    }
    df <- do.call(SingleR,
                  c(list(test = sce,
                         ref = refs,
                         labels = refs$label,
                         BPPARAM = BPPARAM),
                    SingleRParams))
    
    # return results
    list(cells = topL, refs = refs, labels = df)
}
