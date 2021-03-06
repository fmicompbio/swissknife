#' @title Calculate phasograms (same strand alignment distances).
#'
#' @description Calculate the frequencies of same strand alignment distances,
#'   for example from MNase-seq data to estimate nucleosome repeat length.
#'   Distance calculations are implemented in C++ (\code{\link{calcAndCountDist}})
#'   for efficiency.
#'
#' @author Michael Stadler
#'
#' @param fname \code{character} vector with one or several bam files. If
#'   multiple files are given, distance counts from all will be summed.
#' @param regions \code{GRanges} object. Only alignments falling into these
#'   regions will be used. If \code{NULL} (the default), all alignments are
#'   used.
#' @param rmdup \code{logical(1)} indicating if duplicates should be removed. If
#'   \code{TRUE} (the default), only one of several alignments starting at the
#'   same coordinate is used.
#' @param dmax \code{numeric(1)} specifying the maximal distance between same
#'   strand alignments to count.
#'
#' @return \code{integer} vector with \code{dmax} elements, with the element at
#'   position \code{d} giving the observed number of alignment pairs at that
#'   distance.
#'
#' @references Phasograms were originally described in Valouev et al., Nature
#'   2011 (doi:10.1038/nature10002). The implementation here differs in two ways
#'   from the original algorithms:
#'   \enumerate{
#'     \item It does not implement removing of positions that have been seen less
#'     than \code{n} times (referred to as a \code{n}-pile subset in the paper).
#'     \item It does allow to retain only alignments that fall into selected
#'     genomic intervals (\code{regions} argument).
#'   }
#'
#' @seealso \code{\link{estimateNRL}} to estimate the nucleosome repeat length
#'   from a phasogram, \code{\link{plotPhasogram}} to visualize an annotated
#'   phasogram, \code{\link{calcAndCountDist}} for low-level distance counting.
#'
#' @examples
#' if (requireNamespace("GenomicAlignments", quietly = TRUE) &&
#'     requireNamespace("Rsamtools", quietly = TRUE)) {
#'     bamf <- system.file("extdata", "phasograms", "mnase_mm10.bam",
#'                         package = "swissknife")
#'     pg <- calcPhasogram(bamf)
#'     print(estimateNRL(pg, usePeaks = 1:4)[1:2])
#'     plotPhasogram(pg, usePeaks = 1:4, xlim = c(0,1000))
#' }
#'
#' @useDynLib swissknife
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom GenomicRanges GRanges seqnames ranges
#' @import Rcpp
#'
#' @export
calcPhasogram <- function(fname, regions=NULL, rmdup=TRUE, dmax=3000L) {
    .assertPackagesAvailable(c("GenomicAlignments", "Rsamtools"))

    cnt <- numeric(dmax)
    names(cnt) <- as.character(seq.int(dmax))

    # get chromosomes from bam header
    bh <- Rsamtools::scanBamHeader(fname[1])[[1]]$targets
    chrs <- names(bh)
    chrsReg <- GenomicRanges::GRanges(names(bh), IRanges::IRanges(start = 1L, width = bh, names = names(bh)))
    if (!is.null(regions)) {
        stopifnot(inherits(regions, "GRanges"))
        regL <- split(regions, GenomicRanges::seqnames(regions))
        chrs <- intersect(chrs, names(regL))
        regL <- regL[chrs]
    }

    # for each chromosome, ...
    for (i in seq_along(chrs)) {
        posPL <- posML <- list()
        for (j in seq_along(fname)) { # ... and each bam file ...
            # ... read alignment positions
            posP <- Rsamtools::scanBam(fname[j],
                                       param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = FALSE),
                                                                       what = c("pos"), which = chrsReg[chrs[i]]))[[1]]$pos
            alnM <- Rsamtools::scanBam(fname[j],
                                       param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand = TRUE),
                                                                       what = c("pos","cigar"), which = chrsReg[chrs[i]]))[[1]]
            posM <- alnM$pos + GenomicAlignments::cigarWidthAlongReferenceSpace(alnM$cigar, flag = NULL, N.regions.removed = FALSE) - 1L
            posPL <- c(posPL, list(posP))
            posML <- c(posML, list(posM))
        }
        # ... combine
        posP <- do.call(c, posPL)
        posM <- do.call(c, posML)
        # ... make unique (if rmdup)
        if (rmdup) {
            posP <- unique(posP)
            posM <- unique(posM)
        }
        # ... sort
        posP <- sort(posP)
        posM <- sort(posM)
        # ... filter by selected regions
        if (!is.null(regions)) {
            posP <- posP[ IRanges::overlapsAny(IRanges::IRanges(posP, width = 1L), GenomicRanges::ranges(regL[[i]])) ]
            posM <- posM[ IRanges::overlapsAny(IRanges::IRanges(posM, width = 1L), GenomicRanges::ranges(regL[[i]])) ]
        }
        # ... phasogram (alignments on same strand)
        cnt <- calcAndCountDist(posP, posP, cnt) # plus strand
        cnt <- calcAndCountDist(posM, posM, cnt) # minus strand
    }

    return(cnt)
}

#' @title Estimate the nucleosome repeat length (NRL) from a phasogram.
#'
#' @description Estimate the nucleosome repeat length (NRL) from the frequencies
#'   of same-strand alignment distances (phasogram), e.g. generated by
#'   \code{\link{calcPhasogram}}. The NRL is obtained from the slope of a linear
#'   fit to the modes in the phasogram.
#'
#' @author Michael Stadler
#'
#' @param x \code{numeric} vector giving the counts of alignment distances
#'   (typically the output of \code{\link{calcPhasogram}}.
#' @param mind \code{integer(1)} specifying the minimal distance to be used for
#'   NRL estimation. The default value (140) ignores any distance too short to
#'   span at least a single nucleosome.
#' @param usePeaks \code{integer} vector selecting the modes (peaks) in the
#'   phasogram used in NRL estimation.
#' @param span1 \code{numeric(1)} giving the smoothing parameter for de-trending
#'   loess fit (high pass filter).
#' @param span2 \code{numeric(1)} giving the smoothing parameter for de-noising
#'   loess fit (low pass filter).
#'
#' @return A \code{list} with elements: \describe{
#'   \item{nrl}{the estimated nucleosome repeat length}
#'   \item{nrl.CI95}{the 95\% confidence interval}
#'   \item{xs}{smoothed (de-trended) phasogram}
#'   \item{loessfit}{the de-noising fit to the de-trended phasogram}
#'   \item{lmfit}{the linear fit to the phasogram peaks}
#'   \item{peaks}{the peak locations}
#'   \item{mind}{minimal distance included in the fit}
#'   \item{span1}{smoothing parameter for de-trending loess fit}
#'   \item{span2}{smoothing parameter for de-noising loess fit}
#'   \item{usePeaks}{the peaks used in the fit}}
#'
#' @seealso \code{\link{calcPhasogram}} to calculate the phasogram from
#'   alignments, \code{\link{plotPhasogram}} to visualize an annotated phasogram
#'
#' @examples
#'   # see the help for calcPhasogram() for a full example
#'
#' @importFrom stats loess lm confint residuals predict coefficients
#' @importFrom IRanges IRanges Views viewApply
#' @importFrom methods as
#'
#' @export
estimateNRL <- function(x, mind = 140L, usePeaks=1:8, span1=100/length(x),
                        span2=1500/length(x)) {
    stopifnot(is.numeric(x))

    if (all(x == 0)) {
        warning("NRL not estimated (phasogram contains only zeros)")
        return(list(nrl = NA, nrl.CI95 = NA, xs = NA, loessfit = NA, lmfit = NA,
                    peaks = NA, mind = mind, span1 = span1, span2 = span2,
                    usePeaks = usePeaks))
    }

    pos <- seq_along(x)
    xs <- predict(loess(x ~ pos, subset = pos > mind, span = span1), pos)
    fit <- loess(xs ~ pos, subset = pos > mind, span = span2)
    rx <- residuals(fit)
    irpos <- as(rx >= 0, "IRanges")
    xposmax <- IRanges::viewApply(IRanges::Views(rx, irpos), function(y) which.max(as.vector(y))) + mind + start(irpos) - 1
    if (any(!usePeaks %in% seq_along(xposmax))) {
        warning("less peaks detected than selected by `usePeaks`")
        usePeaks <- intersect(usePeaks, seq_along(xposmax))
    }
    lmfit <- stats::lm(xposmax ~ seq_along(xposmax), subset = usePeaks)
    suppressWarnings(cilmfit <- stats::confint(lmfit)) # will warn if a single peak and slope = NA
    list(nrl = unname(coefficients(lmfit)[2]), nrl.CI95 = cilmfit[2,],
         xs = xs, loessfit = fit, lmfit = lmfit,
         peaks = xposmax, mind = mind, span1 = span1, span2 = span2, usePeaks = usePeaks)
}


#' @title Plot annotated phasogram.
#'
#' @description Plot phasogram and annotate it with estimated nucleosome repeat
#'   length (NRL).
#'
#' @author Michael Stadler
#'
#' @param x \code{numeric} vector giving the counts of alignment distances
#'   (typically the output of \code{\link{calcPhasogram}}.
#' @param hide If \code{TRUE} (the default), hide phasogram counts not used in
#'   the NRL estimate (\code{mind} parameter from \code{\link{estimateNRL}}).
#' @param xlim \code{numeric(2)} with the x-axis (phase) limits in the first two
#'   plots (see Details). if \code{NULL} (the default), the full range defined
#'   by \code{x} and \code{hide} will be used.
#' @param verbosePlot If \code{TRUE}, create three plots instead of just a
#'   single plot (see Details).
#' @param ... Additional arguments passed to \code{\link{estimateNRL}} to
#'   control NRL estimation.
#'
#' @details The function will visualize an annotated phasogram. For
#'   \code{verbosePlot=FALSE} (the default), it will create a single annotated
#'   plot. For \code{verbosePlot=TRUE}, it will create three plots (using
#'   \code{par(mfrow=c(1,3))}): \enumerate{ \item raw phase counts with
#'   de-trending and de-noising loess fits \item residual phase counts with
#'   de-noising loess fit and detected peaks \item linear fit to peaks and NRL
#'   estimation }
#'
#' @return The return value from the call to \code{\link{estimateNRL}}
#'   (invisibly).
#'
#' @seealso \code{\link{calcPhasogram}} to calculate the phasogram from
#'   alignments, \code{\link{estimateNRL}} to estimate nucleosome repeat length
#'
#' @examples
#'   # see the help for calcPhasogram() for a full example
#'
#' @importFrom stats residuals
#' @importFrom IRanges IRanges start end
#' @importFrom methods as
#' @import graphics
#'
#' @export
plotPhasogram <- function(x, hide = TRUE, xlim = NULL, verbosePlot = FALSE, ...) {
    nrl <- estimateNRL(x, ...)
    pos <- seq_along(x)
    if (is.null(xlim))
        xlim <- c(0, length(x))
    if (hide) {
        x[seq.int(nrl$mind - 1)] <- NA
        xlim[1] <- nrl$mind
    }

    if (verbosePlot) {
        par(mfrow = c(1, 3), mar = c(5, 4, 4 - 3, 2) + 0.1, cex = 1.25)
        # raw phase counts
        matplot(pos, cbind(x, nrl$xs, c(rep(NA, nrl$mind), nrl$loessfit$fitted)),
                type = "l", lty = 1, lwd = c(1,2,2), xlim = xlim,
                xlab = "Phase (bp)", ylab = "Phase counts", col = c("gray", "red", "green3"), axes = FALSE)
        axis(1); axis(2)
        pusr <- par('usr')
        segments(x0 = pusr[c(1,1)], y0 = pusr[c(4,3)], x1 = pusr[c(1,2)], y1 = pusr[c(3,3)])
        legend("topright", c("raw", paste0("smoothed (",signif(c(nrl$span1, nrl$span2),4),")")),
               lty = 1, lwd = c(1,2,2), col = c("gray", "red", "green3"), bty = "n")
        # residual phase counts
        rx <- stats::residuals(nrl$loessfit)
        irpos <- as(rx >= 0, "IRanges")
        plot(pos, c(rep(NA, nrl$mind), rx), type = "l", axes = FALSE,
             xlim = xlim, xlab = "Phase (bp)", ylab = "Residual phase counts")
        abline(h = 0, lty = 2)
        axis(1); axis(2)
        pusr <- par('usr')
        segments(x0 = pusr[c(1,1)], y0 = pusr[c(4,3)], x1 = pusr[c(1,2)], y1 = pusr[c(3,3)])
        rect(nrl$mind + IRanges::start(irpos), pusr[3], nrl$mind + IRanges::end(irpos),
             pusr[4], col = "#FF000022", border = NA)
        text(x = nrl$mind + IRanges::start(irpos)[1] - par("cxy")[1],
             y = pusr[4] + 0.25 * par("cxy")[2], adj = c(1,0), labels = "peaks:", xpd = NA)
        points(nrl$peaks, c(rep(NA, nrl$mind), rx)[nrl$peaks], pch = 20)
        text(x = nrl$peaks[nrl$usePeaks], pusr[4] + 0.25 * par("cxy")[2], adj = c(0.5,0),
             labels = nrl$usePeaks, xpd = NA)
        # linear fit
        slmfit <- stats::summary.lm(nrl$lmfit)
        plot(seq_along(nrl$peaks)[nrl$usePeaks], nrl$peaks[nrl$usePeaks], pch = 21,
             xlab = "Peak count", ylab = "Peak coordinate (bp)")
        abline(nrl$lmfit)
        legend(x = "topleft", bty = "n",
               legend = sprintf("%1g (%1g-%1g)", signif(nrl$nrl,3),
                                signif(nrl$nrl.CI95[1],3), signif(nrl$nrl.CI95[2],3)))
        legend(x = "bottomright", bty = "n",
               legend = c(paste0("Adj. R-squared: ", signif(slmfit$adj.r.squared, 3)),
                          paste0("P-value = ", signif(slmfit$coefficients[2,4],3))))
    } else {
        par(mar = c(5, 4, 4 - 3, 2) + 0.1, cex = 1.25)
        # raw phase counts
        matplot(pos, cbind(x, nrl$xs, c(rep(NA, nrl$mind), nrl$loessfit$fitted)),
                type = "l", lty = 1, lwd = c(1,2,2), xlim = xlim,
                xlab = "Phase (bp)", ylab = "Phase counts", col = c("gray", "red", "green3"), axes = FALSE)
        axis(1); axis(2)
        pusr <- par('usr')
        segments(x0 = pusr[c(1,1)], y0 = pusr[c(4,3)], x1 = pusr[c(1,2)], y1 = pusr[c(3,3)])
        legend("bottomleft", c("raw", paste0("smoothed (",signif(c(nrl$span1, nrl$span2),4),")")),
               lty = 1, lwd = c(1,2,2), col = c("gray", "red", "green3"), bty = "n")
        rx <- stats::residuals(nrl$loessfit)
        irpos <- as(rx >= 0, "IRanges")
        text(x = nrl$mind + IRanges::start(irpos)[1] - par("cxy")[1],
             y = pusr[4] + 0.25 * par("cxy")[2], adj = c(1,0), labels = "peaks:", xpd = NA)
        points(nrl$peaks[nrl$usePeaks], nrl$xs[nrl$peaks][nrl$usePeaks], pch = 20, col = "red")
        text(x = nrl$peaks[nrl$usePeaks], pusr[4] + 0.25 * par("cxy")[2], adj = c(0.5,0),
             labels = nrl$usePeaks, xpd = NA)
        slmfit <- stats::summary.lm(nrl$lmfit)
        legend(x = "topright", bty = "n",
               legend = sprintf("%1g (%1g-%1g)", signif(nrl$nrl,3),
                                signif(nrl$nrl.CI95[1],3), signif(nrl$nrl.CI95[2],3)))
    }
    return(invisible(nrl))
}

