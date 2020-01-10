#' @title Get Distance
#' 
#' @description This function gets the euclidean (by default) distances between all points in a 2-D coordinate 
#'   system (x,y) to a group of points that fall on a fitted loess curve. It uses the \code{\link[wordspace]{dist.matrix}}
#'   function from the \code{wordspace} package.It is used in the \code{selVarGenes} function.
#' 
#' @author Dania Machlab 
#' 
#' @param loessModel the loess fit to all points.
#' @param x the x values of all points.
#' @param y the y values of all points.
#' @param x_curve the x values of the points on the loess curve.
#' @param method distance measure to be used.
#' @param ... additional parameters for the \code{dist.matrix} function.
#' 
#' @return An object of class \code{matrix} and \code{dist.matrix}, containing the
#'    distances between all points to the specified points on the loess curve.
#'
#' @examples 
#'   library(stats)
#'   
#'   # create data
#'   x <- 1:1000
#'   y <- sort(rnorm(n = 1000, mean = 0, sd = 5)) + runif(n = 1000, min = 0, max = 10)
#'   lo <- loess(y ~ x)
#'   x_curve <- seq(range(x)[1], range(x)[2], length.out = 100)
#'   
#'   # plot
#'   plot(x, y, pch = 16, col = "#BEBEBE90")
#'   lines(x_curve, predict(lo, newdata = x_curve), col = "red")
#'   
#'   # get euclidean distance of all points to the curve (red)
#'   euclMat <- swissknife:::.getDistMat(loessModel = lo, x = x, y = y, x_curve = x_curve)
#'   euclMat[1:6, 1:6]
#'
#' @importFrom wordspace dist.matrix
#' @importFrom stats predict
#' 
#' @keywords internal
#' 
.getDistMat <- function(loessModel = NULL, x = NULL, y = NULL, x_curve = NULL, method = "euclidean", ...){
     
     ## checks
     stopifnot(!is.null(x_curve))
     stopifnot(!is.null(x))
     stopifnot(!is.null(y))
     stopifnot(!is.null(loessModel))
     stopifnot(length(x) == length(y))
     stopifnot(is.numeric(x_curve))
     stopifnot(is.numeric(x))
     stopifnot(is.numeric(y))

     ## predict y_curve for x_curve using loessModel
     y_curve <- stats::predict(loessModel, newdata = x_curve)
     
     ## get distances
     distMat <- wordspace::dist.matrix(M = cbind(x, y), M2 = cbind(x_curve, y_curve), method = method, ... )
     colnames(distMat) <- paste0("pointOnCurve", seq_len(length(y_curve)))
     if (!is.null(names(x))) {
          rownames(distMat) <- names(x)
     } else {
          rownames(distMat) <- paste0("xy", seq_len(length(x)))
     }
     
     ## return
     distMat
     
}


#' @title Select Variable Genes in Single Cell RNA-seq
#' 
#' @description This function selects the most variable genes from a \code{SingleCellExperiment} object
#'   using the plot that displays the log2 coefficient of variation as a function of the log2 mean for 
#'   all genes across all the cells.
#' 
#' @author Dania Machlab 
#' 
#' @param data \code{SingleCellExperiment} object or (normalized) count \code{matrix} containing the genes as rows and cells as columns.
#' @param assay.type the type of assay to use if \code{data} is a \code{SingleCellExperiment}. It can be either 'counts' 
#'    or 'logcounts'. The default is 'counts'.
#' @param Nmads number of MADs beyond which genes are selected per bin.
#' @param minCells keep genes with minimum expression in at least this number of cells.
#' @param minExpr keep genes with expression greater than or equal to this in \code{minCells} cells in the 
#'    normalized count matrix.
#' @param exclTopExprFrac the fraction of top expressed genes that will be excluded from the loess fit 
#'    (value between 0 and 1). 
#' @param span span parameter for \code{loess} function.
#' @param control control parameters for \code{loess} function.
#' @param nBins number of bins or groups to place the points(genes) into.
#' @param nBinsDense number of bins or groups to use  to place the points(genes) into
#'    when calculating more accurate distance values to the curve from the loess fit.
#' @param ... additional parameters for the \code{loess} function from the \code{stats} package.
#'
#' @details The function takes in a \code{SingleCellExperiment} object and calculates the normalized
#'   counts by dividing the raw counts by the corresponding sizeFactors per cell. Only genes that have
#'   an expression greater than or equal to \code{minExpr} in at least \code{minCells} cells will be
#'   kept.
#'   
#'   The genes that vary most on the log2(coefficient of variation) vs log2(mean) plot of genes will 
#'   be selected. A loess fit is done on this plot and the distance (euclidean by default) each point has
#'   to the curve is calculated in two steps.
#'   
#'   In the first step, genes are assigned to bins by taking the minimum distance to the curve. By default
#'   we select 100 points on the loess fit and calculate the distances each gene has to all those points on 
#'   the curve. Each gene is assigned to the point on the curve for which it has the shortest distance. In 
#'   the second step, more accurate distances to the curve are calculated by using a higher number of 
#'   points on the curve. Distances are calculated using the \code{\link[wordspace]{dist.matrix}} function.
#'   
#'   Finally, for each bin, the most variable genes are selected using the more accurate distance
#'   measures. Genes that fall below the loess fit are assigned a negative sign and the genes that 
#'   are \code{Nmads} MADs away from the median are selected.
#'   
#' 
#' @return a \code{list} of length 2:
#'   \itemize{
#'      \item varGenes: \code{vector} containing the names of the most variable genes.
#'      \item geneInfo: \code{data.frame} with genes as rows and columns containing  
#'                      calculated measures:
#'                      \itemize{
#'                         \item logMean: log2(mean) expression of genes across cells.
#'                         \item logCV: log2(coefficient of variation) of genes across cells.
#'                         \item pred_logCV: predicted log2(coefficient of variation) from loess fit.
#'                         \item assigned_bin: bin each gene has been assigned to.
#'                         \item distance: second round of calculating more accurate distances
#'                      }
#'   }
#'   
#' @examples 
#'    # packages
#'    library(SingleCellExperiment)
#'    
#'    # create example count matrix
#'    # ... poisson distr per gene
#'    mu <- ceiling(runif(n = 2000, min = 0, max = 100))
#'    counts <- do.call(rbind, lapply(mu, function(x){rpois(1000, lambda = x)}))
#'    counts <- counts + 1
#'    # ... add signal to subset of genes (rows) and cells (columns)
#'    i <- sample(x = 1:nrow(counts), size = 500)
#'    j <- sample(x = 1:ncol(counts), size = 500)
#'    counts[i, j] <- counts[i, j] + sample(5:10, length(i), replace = TRUE)
#'    
#'    # create SCE 
#'    sce <- SingleCellExperiment(list(counts=counts))
#'    
#'    # calculate sizeFactors
#'    libsizes <- colSums(counts)
#'    sizeFactors(sce) <- libsizes/mean(libsizes)
#'  
#'    # select variable genes
#'    varGenes <- selVarGenes(sce, assay.type="counts")
#'    
#'    # plot
#'    plotSelVarGenes(varGenes, colByBin=TRUE)
#'    plotSelVarGenes(varGenes)
#'
#' @importFrom stats loess loess.control predict mad median quantile sd
#' @importFrom SingleCellExperiment sizeFactors counts logcounts
#' @importFrom SummarizedExperiment assays
#' 
#' @export
#' 
selVarGenes <- function(data=NULL, assay.type="counts", Nmads = 3, minCells = 5, minExpr = 1, exclTopExprFrac = 0.01, span = 0.2, 
                        control=stats::loess.control(surface = "direct"), nBins = 100, 
                        nBinsDense = ceiling(nrow(data)/4), ...){
     
     ## checks
     if (is.null(data)) {stop("'data' is empty")}
     if (!is(data, "SingleCellExperiment") & !is(data, "matrix")) {stop("data must be a SingleCellExperiment object or a count matrix")}
     if (is(data, "SingleCellExperiment") & !any(assay.type%in%c("counts", "logcounts"))) {stop("assay.type must be either counts or logcounts")}
     if (is(data, "SingleCellExperiment") & assay.type=="counts") {
          if(!any(names(SummarizedExperiment::assays(data)) %in% "counts")) {stop("counts assay is missing from SCE object")}
     }
     if (is(data, "SingleCellExperiment") & assay.type=="logcounts") {
          if(!any(names(SummarizedExperiment::assays(data)) %in% "logcounts")) {stop("logcounts assay is missing from SCE object")}
     }
     if (is(data, "SingleCellExperiment") & assay.type=="counts") {
          if(is.null(SingleCellExperiment::sizeFactors(data))) {stop("missing sizeFactors to normalize counts")}
     }
        
     ## when data is SCE
     if(is(data, "SingleCellExperiment")) {
          if(assay.type=="logcounts"){
               normCounts <- SingleCellExperiment::logcounts(data)
          } else {
               normCounts <- sweep(as.matrix(SingleCellExperiment::counts(data)), 2, SingleCellExperiment::sizeFactors(data), "/")  
          }
             
     } 
     
     ## when data is count matrix
     if(is(data, "matrix")){
          message("assuming that data is a normalized count matrix with genes as rows and cells as columns, if that's not the case make sure it is ... ")
          normCounts <- data
     }
        
     ## normCounts rownames and colnames
     if (is.null(rownames(normCounts))) {
          rownames(normCounts) <- paste0("Gene_", seq_len(nrow(normCounts)))
          message("Count row names are empty, naming them now ... ")
     }
     if (is.null(colnames(normCounts))) {
          colnames(normCounts) <- paste0("Cell_", seq_len(ncol(normCounts)))
          message("Count column names are empty, naming them now ... ")
     }
     
     ## keep genes with minExpression in at least minCells
     keep <- rowSums(normCounts >= minExpr) >= minCells
     normCounts <- normCounts[keep, , drop = FALSE]
     
     ## log2(mean expression) and log2(coefficient of variation)
     logMean <- log2(rowMeans(normCounts))
     logCV <- log2(apply(normCounts, 1, stats::sd)/rowMeans(normCounts))
     datfr <- data.frame(logMean = logMean, logCV = logCV)
     
     ## loess fit, excluding the top exclTopExprFrac
     points_for_loess <- which(logMean < stats::quantile(logMean, (1 - exclTopExprFrac)))
     lo <- stats::loess(logCV ~ logMean, span = span, control = control, 
                        subset = points_for_loess, ...)
     datfr$pred_logCV <- stats::predict(lo, newdata = datfr$logMean)
     
     ## assign to groups (min eucl distance to points on the loess curve)
     dist1 <- .getDistMat(loessModel = lo, x = logMean, y = logCV, 
                          x_curve = seq(from = range(logMean)[1], 
                                        to = range(logMean)[2], 
                                        length.out = nBins))
     datfr$assigned_bin <- apply(X = dist1, MARGIN = 1, FUN = function(x){which.min(x)}) 
     # datfr$min_dist1_per_gene <- sapply(seq_along(datfr$assigned_bin), 
     #                                    function(i){dist1[i, datfr$assigned_bin[i]]})
      
     ## get more accurate eucl dist to the curve for each gene
     dist2 <- .getDistMat(loessModel = lo, x = logMean, y = logCV, 
                         x_curve = seq(from = range(logMean)[1], 
                                       to = range(logMean)[2], 
                                       length.out = nBinsDense))
     sel <- apply(X = dist2, MARGIN = 1, FUN = function(x){which.min(x)}) 
     datfr$distance <- sapply(seq_along(sel), function(i){dist2[i, sel[i]]})
     
     ## give dist2 a sign: + for above the curve, and - for below the curve
     y_pred_for_x <- stats::predict(lo, newdata = logMean)
     w_down <- logCV < y_pred_for_x
     datfr$distance[w_down] <- -datfr$distance[w_down]
     
     ## genes per group
     genes_per_bin <- lapply(unique(datfr$assigned_bin), 
                             function(i){rownames(datfr)[datfr$assigned_bin == i]})
     names(genes_per_bin) <- unique(datfr$assigned_bin)
     
     ## MAD (this includes genes that were excluded from the loess) per group
     mads <- sapply(genes_per_bin, function(x){stats::mad(datfr[x, ]$distance)})
     mads <- mads*Nmads
     
     ## select outlier genes per group
     out_genes <- unlist(lapply(seq_along(genes_per_bin), function(i){
          genes <- genes_per_bin[[i]]
          distances <- datfr[genes, ]$distance
          median <- stats::median(distances)
          genes[distances > (median + mads[i])]
     }))
     out_genes <- out_genes[!is.na(out_genes)]
     
     ## return outlier genes
     list(varGenes = out_genes, geneInfo = datfr)
     
}


#' @title Plot Selected Variable Genes
#' 
#' @description This function take the output from \code{selVarGenes} and plots the genes that have been 
#'   selected to be highly variable across the cells. It plot the log2 coefficient of variation as
#'   a function of the log mean.
#' 
#' @author Dania Machlab 
#' 
#' @param selVarGenes_list the output list from the \code{selVarGenes} function.
#' @param xlab label for x-axis.
#' @param ylab label for y-axis.
#' @param main title for plot.
#' @param pch point pch.
#' @param col point color.
#' @param sel_col point color of the selected variable genes.
#' @param colByBin if TRUE, color the genes by the bin they've been assigned to.
#' @param asp the y/x aspect ratio. Set to 1 when \code{colByBin} is TRUE.
#' @param ... additional parameters for the \code{plot} function.
#' 
#' @return plot
#'
#' @examples 
#'    # packages
#'    library(SingleCellExperiment)
#'    
#'    # create example count matrix
#'    # ... poisson distr per gene
#'    mu <- ceiling(runif(n = 2000, min = 0, max = 100))
#'    counts <- do.call(rbind, lapply(mu, function(x){rpois(1000, lambda = x)}))
#'    counts <- counts + 1
#'    # ... add signal to subset of genes (rows) and cells (columns)
#'    i <- sample(x = 1:nrow(counts), size = 500)
#'    j <- sample(x = 1:ncol(counts), size = 500)
#'    counts[i, j] <- counts[i, j] + sample(5:10, length(i), replace = TRUE)
#'    
#'    # create SCE
#'    sce <- SingleCellExperiment(list(counts=counts))
#'    
#'    # calculate sizeFactors
#'    libsizes <- colSums(counts)
#'    sizeFactors(sce) <- libsizes/mean(libsizes)
#'  
#'    # select variable genes
#'    varGenes <- selVarGenes(sce)
#'    
#'    # plot
#'    plotSelVarGenes(varGenes)
#'    plotSelVarGenes(varGenes, colByBin=TRUE)
#' 
#' @importFrom grDevices colors
#' 
#' @export
#' 
plotSelVarGenes <- function(selVarGenes_list = NULL, xlab = "logMean", 
                            ylab = "logCV", main = "Selected Variable Genes", 
                            pch = 16, col = "#BEBEBE40", sel_col = "steelblue", 
                            colByBin = FALSE, asp = 1, ...) {
     
     ## checks
     if (any(is.null(selVarGenes_list))) {
          stop("'selVarGenes_list' is empty")
     }
     if (!is(selVarGenes_list, "list")) {
          stop("'selVarGenes_list' must be of class 'list'")
     }
     if (!all(names(selVarGenes_list) %in% c("varGenes", "geneInfo"))) {
          stop("names of 'selVarGenes_list' must be 'varGenes' and 'geneInfo', 
           the output from the 'selVarGenes' function")
     }
             
     ## plot
     if(colByBin) {
             
          ## prepare for plot
          colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
          yHat_cols <- sample(colors, length(unique(selVarGenes_list$geneInfo$assigned_bin)))
          grpCols <- yHat_cols[selVarGenes_list$geneInfo$assigned_bin]
             
          ## plot
          plot(selVarGenes_list$geneInfo$logMean, selVarGenes_list$geneInfo$logCV, 
               xlab = xlab, ylab = ylab, main = main, pch = pch, col = grpCols, 
               asp = asp, ...)
          
          ## return TRUE
          invisible(TRUE)
             
     } else {
        
          plot(selVarGenes_list$geneInfo$logMean, selVarGenes_list$geneInfo$logCV, 
               xlab = xlab, ylab = ylab, main = main, pch = pch, col = col, ...)
          points(selVarGenes_list$geneInfo[selVarGenes_list$varGenes, ]$logMean, 
                 selVarGenes_list$geneInfo[selVarGenes_list$varGenes, ]$logCV, 
                 pch = pch, col = sel_col)
          lines(x = selVarGenes_list$geneInfo$logMean[order(selVarGenes_list$geneInfo$logMean)], 
                y = selVarGenes_list$geneInfo$pred_logCV[order(selVarGenes_list$geneInfo$logMean)], 
                col = "red")
          legend("bottomleft", bty = "n", col = c(sel_col, "red"), pch = c(pch, NA), 
                 lty = c(NA, 1), legend = c("SelVarGene", "loess fit"))
          
          ## return TRUE
          invisible(TRUE)
     
     }

     
}



