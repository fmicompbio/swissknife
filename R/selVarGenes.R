#' @title Get Distance
#' 
#' @description This function gets the euclidean (by default) distances between all points in a 2-D coordinate 
#'   system (x,y) to a group of points that fall on a fitted loess curve. It uses the \code{dist.matrix} 
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
#'   x_curve <- seq(range(x)[1], range(x)[2], length.out=100)
#'   
#'   # plot
#'   plot(x, y, pch = 16, col="#BEBEBE90")
#'   lines(x_curve, predict(lo, newdata=x_curve), col="red")
#'   
#'   # get euclidean distance of all points to the curve (red)
#'   euclMat <- getDistMat(loessModel=lo, x=x, y=y, x_curve=x_curve)
#'   euclMat[1:6, 1:6]
#'
#' @importFrom wordspace dist.matrix
#' 
#' @export
#' 
getDistMat <- function(loessModel=NULL, x=NULL, y=NULL, x_curve=NULL, method="euclidean", ...){
     
     ## checks
     stopifnot(!is.null(x_curve))
     stopifnot(!is.null(x))
     stopifnot(!is.null(y))
     stopifnot(!is.null(loessModel))
     stopifnot(length(x)==length(y))
     stopifnot(is.numeric(x_curve))
     stopifnot(is.numeric(x))
     stopifnot(is.numeric(y))
     
     ## predict y_curv for x_curve using loessModel
     y_curve <- predict(loessModel, newdata = x_curve)
     
     ## get distances
     distMat <- wordspace::dist.matrix(M = cbind(x,y), M2 = cbind(x_curve, y_curve), method = method, ... )
     colnames(distMat) <- paste0("pointOnCurve", 1:length(y_curve))
     if(!is.null(names(x))){
          rownames(distMat) <- names(x)
     } else {
          rownames(distMat) <- paste0("xy", 1:length(x))
     }
     
     ## return
     distMat
     
}


#' @title Select Variable Genes in Single Cell RNA-seq
#' 
#' @description This function selects the most variable genes from a \code{Single Cell Experiment} object
#'   using the plot that displays the log2 coefficient of variation as a function of the log2 mean for 
#'   all genes across all the cells.
#' 
#' @author Dania Machlab 
#' 
#' @param sce \code{SingleCellExperiment} object containing the raw counts and the scaling factors.
#' @param Nmads number of MADs beyond which genes are selected per bin.
#' @param minCells keep genes with minimum expression in at least this number of cells.
#' @param minExpr keep genes with expression greater than or equal to this in \code{minCells} cells using 
#'    the raw count matrix corrected for library size differences using the sizeFactors from \code{sce}.
#' @param topExprPerc the top percent expressed genes that will be excluded from the loess fit. 
#' @param span span parameter for \code{loess} function.
#' @param control control parameters for \code{loess} function.
#' @param binBreaks number of bins or groups to place the points(genes) into.
#' @param accurateDistBreaks number of bins or groups to use  to place the points(genes) into
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
#'   points on the curve.
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
#'                         \item bin_per_gene: bin each gene has been assigned to.
#'                         \item min_dist1_per_gene: first round of calcultaing distances.
#'                         \item min_dist2_per_gene: second round of calcultaing more accurate distances
#'                      }
#'   }
#'   
#' @examples 
#'
#' @importFrom stats loess loess.control predict mad median quantile sd
#' @importFrom SingleCellExperiment sizeFactors counts
#' 
#' @export
#' 
selVarGenes <- function(sce=NULL, Nmads = 3, minCells = 5, minExpr = 1, topExprPerc = 0.01, span = 0.2, 
                        control=stats::loess.control(surface = "direct"), binBreaks=100, 
                        accurateDistBreaks = ceiling(nrow(sce)/4), ...){
     
     ## checks
     if(is.null(sce)){stop("'sce' is empty")}
     if(class(sce)!="SingleCellExperiment"){stop("sce must be a SingleCellExperiment object")}
     if(is.null(SingleCellExperiment::sizeFactors(sce))){stop("sce must contain sizeFactors to normalize the raw counts")}
     
     ## normCounts
     normCounts <- sweep(as.matrix(SingleCellExperiment::counts(sce)), 2, SingleCellExperiment::sizeFactors(sce), "/")
     if(is.null(rownames(normCounts))){
          rownames(normCounts) <- paste0("Gene_", 1:nrow(normCounts))
          message("Count row names are empty, naming them now ... ")
     }
     if(is.null(colnames(normCounts))){
          colnames(normCounts) <- paste0("Cell_", 1:ncol(normCounts))
          message("Count column names are empty, naming them now ... ")
     }
     
     ## keep genes with minExpression in at least minCells
     keep <- rowSums(normCounts>=minExpr) >= minCells
     normCounts <- normCounts[keep, ]
     
     ## log2(mean expression) and log2(coefficient of variaition)
     logMean <- log2(rowMeans(normCounts))
     logCV <- log2(apply(normCounts, 1, stats::sd)/rowMeans(normCounts))
     datfr <- data.frame(logMean = logMean, logCV = logCV)
     
     ## loess fit, excluding the top topExprPerc
     points_for_loess <- which(logMean < stats::quantile(logMean, (1-topExprPerc)))
     lo <- stats::loess(logCV ~ logMean, span = span, control = control, subset = points_for_loess, ...)
     datfr$pred_logCV <- stats::predict(lo, newdata = datfr$logMean)
     
     ## assign to groups (min eucl distance to points on the loess curve)
     dist1 <- getDistMat(loessModel=lo, x=logMean, y=logCV, x_curve=seq(from = range(logMean)[1], to = range(logMean)[2], length.out = binBreaks))
     datfr$bin_per_gene <- apply(X = dist1, MARGIN = 1, FUN = function(x){which.min(x)}) 
     datfr$min_dist1_per_gene <- sapply(seq_along(datfr$bin_per_gene), function(i){dist1[i, datfr$bin_per_gene[i]]})
     
     ## get more accurate eucl dist to the curve for each gene
     dist2 <- getDistMat(loessModel=lo, x=logMean, y=logCV, x_curve=seq(from = range(logMean)[1], to = range(logMean)[2], length.out = accurateDistBreaks))
     sel <- apply(X = dist2, MARGIN = 1, FUN = function(x){which.min(x)}) 
     datfr$min_dist2_per_gene <- sapply(seq_along(sel), function(i){dist2[i, sel[i]]})
     
     ## give dist2 a sign: + for above the curve, and - for below the curve
     y_pred_for_x <- stats::predict(lo, newdata = logMean)
     w_down <- logCV < y_pred_for_x
     datfr$min_dist2_per_gene[w_down] <- -datfr$min_dist2_per_gene[w_down]
     
     ## genes per group
     genes_per_bin <- lapply(unique(datfr$bin_per_gene), function(i){rownames(datfr)[datfr$bin_per_gene==i]})
     names(genes_per_bin) <- unique(datfr$bin_per_gene)
     
     ## MAD (this includes genes that were excluded from the loess) per group
     mads <- sapply(genes_per_bin, function(x){stats::mad(datfr[x, ]$min_dist2_per_gene)})
     mads <- mads*Nmads
     
     ## select outlier genes per group
     out_genes <- unlist(lapply(seq_along(genes_per_bin), function(i){
          genes <- genes_per_bin[[i]]
          distances <- datfr[genes, ]$min_dist2_per_gene
          median <- stats::median(distances)
          genes[distances > (median + mads[i])]
     }))
     out_genes <- out_genes[!is.na(out_genes)]
     
     ## return outlier genes
     list(varGenes=out_genes, geneInfo=datfr)
     
}



#' @title Plot Groups Genes Have Been Assigned to
#' 
#' @description This function take the output from \code{selVarGenes} and plots the group or bin each
#'   gene has been assigned to on the scatter plot showing the log2 coefficient of variation as a
#'   function of the log2 mean.
#' 
#' @author Dania Machlab 
#' 
#' @param selVarGenes_list the output list from the \code{selVarGenes} function.
#' @param xlab label for x-axis.
#' @param ylab label for y-axis.
#' @param main title for plot.
#' @param asp the y/x aspect ratio.
#' @param ... additional parameters for the \code{plot} function.
#'
#' @details 
#' 
#' @return plot
#'
#' @examples 
#'
#' @importFrom grDevices colors
#' 
#' @export
#' 
plotSelVarGenesGroups <- function(selVarGenes_list = NULL, xlab="logMean", 
                                  ylab="logCV", main = "Group Assignment per Gene", 
                                  asp=1, ...) {
     
     ## checks
     if(any(is.null(selVarGenes_list))){
          stop("'selVarGenes_list' is empty")
     }
     if(class(selVarGenes_list)!="list"){
          stop("'selVarGenes_list' must be of class 'list'")
     }
     if(!all(names(selVarGenes_list)%in%c("varGenes", "geneInfo"))){
          stop("names of 'selVarGenes_list' must be 'varGenes' and 'geneInfo', 
           the output from the 'selVarGenes' function")
     }
     
     ## prepare for plot
     colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert=TRUE)]
     yHat_cols <- sample(colors, length(unique(selVarGenes_list$geneInfo$bin_per_gene)))
     cols <- yHat_cols[selVarGenes_list$geneInfo$bin_per_gene]
     
     ## plot
     plot(selVarGenes_list$geneInfo$logMean, selVarGenes_list$geneInfo$logCV, bg=cols, pch=21, 
          asp=asp, xlab=xlab, ylab=ylab, main=main, ...)
     
     ## return TRUE
     invisible(TRUE)
     
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
#' @param ... additional parameters for the \code{plot} function.
#'
#' @details 
#' 
#' @return plot
#'
#' @examples 
#' 
#' @export
#' 
plotSelVarGenes <- function(selVarGenes_list = NULL, xlab="logMean", 
                            ylab="logCV", main = "Selected Variable Genes", 
                            pch=16, col="#BEBEBE40", sel_col="steelblue", ...) {
     
     ## checks
     if(any(is.null(selVarGenes_list))){
          stop("'selVarGenes_list' is empty")
     }
     if(class(selVarGenes_list)!="list"){
          stop("'selVarGenes_list' must be of class 'list'")
     }
     if(!all(names(selVarGenes_list)%in%c("varGenes", "geneInfo"))){
          stop("names of 'selVarGenes_list' must be 'varGenes' and 'geneInfo', 
           the output from the 'selVarGenes' function")
     }
     
     ## plot
     plot(selVarGenes_list$geneInfo$logMean, selVarGenes_list$geneInfo$logCV, 
          xlab=xlab, ylab=ylab, main=main, pch=pch, col=col, ...)
     points(selVarGenes_list$geneInfo[selVarGenes_list$varGenes, ]$logMean, selVarGenes_list$geneInfo[selVarGenes_list$varGenes, ]$logCV, 
            pch=pch, col=sel_col)
     lines(x = selVarGenes_list$geneInfo$logMean[order(selVarGenes_list$geneInfo$logMean)], 
           y = selVarGenes_list$geneInfo$pred_logCV[order(selVarGenes_list$geneInfo$logMean)], 
           col="red")
     legend("bottomleft", bty="n", col=c(sel_col, "red"), pch=c(pch, NA), 
            lty=c(NA, 1), legend=c("SelVarGene", "loess fit"))
     
     ## return TRUE
     invisible(TRUE)
     
}



