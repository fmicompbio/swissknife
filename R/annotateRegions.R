#' @title Annotate regions.
#'
#' @description Annotate a \code{\link[GenomicRanges]{GRanges}} object with
#'     sets of reference \code{\link[GenomicRanges]{GRanges}} or
#'     \code{\link[GenomicRanges]{GRangesList}} objects, with
#'     respect to overlaps and nearest neighbors.
#'
#' @author Michael Stadler
#'
#' @param x The \code{\link[GenomicRanges]{GRanges}} object to annotate.
#' @param hasOverlap Named \code{list} with \code{\link[GenomicRanges]{GRanges}}
#'     or \code{\link[GenomicRanges]{GRangesList}} object(s).
#'     For each list element, a logical vector "X.hasOverlap" will be added to the
#'     \code{mcols} of the result, with \code{TRUE} for each tile that overlaps
#'     any region in that element. "X" is obtained from \code{names(hasOverlap)}.
#' @param fracOverlap  Named \code{list} with \code{\link[GenomicRanges]{GRanges}}
#'     or \code{\link[GenomicRanges]{GRangesList}} object(s).
#'     For each list element, a numeric vector "X.fracOverlap" will be added to the
#'     \code{mcols} of the result, with a value between 0 and 1 giving the fraction
#'     of bases in a tile that overlaps with any region in that element. "X" is
#'     obtained from \code{names(fracOverlap)}.
#' @param numOverlap  Named \code{list} with \code{\link[GenomicRanges]{GRanges}}
#'     or \code{\link[GenomicRanges]{GRangesList}} object(s).
#'     For each list element, two numeric vectors "X.numOverlapWithin" and
#'     "X.numOverlapAny" will be added to the \code{mcols} of the result, giving
#'     the number of ranges in that element that are fully contained within
#'     a tile, or that overlap with a tile in any way, respectively. "X" is
#'     obtained from \code{names(numOverlap)}.
#' @param nearest Named \code{list} with \code{\link[GenomicRanges]{GRanges}}
#'     or \code{\link[GenomicRanges]{GRangesList}} object(s).
#'     For each list element, two numeric vectors "X.nearestName" and
#'     "X.nearestDistance" will be added to the \code{mcols} of the result, giving
#'     the name and distance of the nearest range in that element for each tile. "X" is
#'     obtained from \code{names(nearest)}, and the values of "X.nearestName" from
#'     \code{names(nearest$X)}. If multiple nearest ranges are at the same
#'     distance from a tile, an arbitrary one is reported in "X.nearestName".
#' @param ignore.strand Logical scalar passed to
#'     \code{\link[GenomicRanges]{findOverlaps}} when searching for overlaps
#'     between \code{x} and reference regions.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} similar to \code{x}, with
#'     annotations added to its metadata columns (\code{mcols}).
#'
#' @examples
#' library(GenomicRanges)
#' 
#' x <- GRanges("chr1", IRanges(c(1, 12), width = 10))
#' tss <- GRanges("chr1", IRanges(c(1, 10, 30), width = 1,
#'                                names = paste0("t", 1:3)))
#' blacklist <- GRanges("chr1", IRanges(20, width = 5))
#' annotateRegions(x, hasOverlap = list(Blacklist = blacklist),
#'                 fracOverlap = list(Blacklist = blacklist),
#'                 numOverlap = list(TSS = tss),
#'                 nearest = list(TSS = tss))
#'
#' @seealso \code{\link{getGenomicTiles}} that uses this function,
#'     \code{\link[GenomicRanges]{findOverlaps}} and
#'     \code{\link[GenomicRanges]{nearest}} in package \pkg{GenomicRanges} used
#'     internally.
#'
#' @importFrom GenomicRanges GRanges width findOverlaps nearest distance
#' @importFrom methods is
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics unlist
#' 
#' @export
annotateRegions <- function(x,
                            hasOverlap = list(),
                            fracOverlap = list(),
                            numOverlap = list(),
                            nearest = list(),
                            ignore.strand = TRUE) {
    ## check arguments
    ## ... genome and addSeqComp
    stopifnot(exprs = {
        # x
        is(x, "GRanges")
        # hasOverlap
        is.list(hasOverlap)
        length(hasOverlap) == 0L || !is.null(names(hasOverlap))
        all(unlist(lapply(hasOverlap, is, "GRanges"))) || all(unlist(lapply(hasOverlap, is, "GRangesList")))
        # fracOverlap
        is.list(fracOverlap)
        length(fracOverlap) == 0L || !is.null(names(fracOverlap))
        all(unlist(lapply(fracOverlap, is, "GRanges"))) || all(unlist(lapply(fracOverlap, is, "GRangesList")))
        # numOverlap
        is.list(numOverlap)
        length(numOverlap) == 0L || !is.null(names(numOverlap))
        all(unlist(lapply(numOverlap, is, "GRanges"))) || all(unlist(lapply(numOverlap, is, "GRangesList")))
        # nearest
        is.list(nearest)
        length(nearest) == 0L || !is.null(names(nearest))
        all(unlist(lapply(nearest, is, "GRanges"))) || all(unlist(lapply(nearest, is, "GRangesList")))
        all(unlist(lapply(nearest, function(x) !is.null(names(x)))))
        # ignore.strand
        is.logical(ignore.strand)
        length(ignore.strand) == 1L
    })
    
    ## annotate tiles
    df <- S4Vectors::mcols(x)
    
    ## ... hasOverlap
    for (i in seq_along(hasOverlap)) {
        nm <- paste0(names(hasOverlap)[i], ".hasOverlap")
        df[[nm]] <- IRanges::overlapsAny(query = x, subject = hasOverlap[[i]],
                                         ignore.strand = ignore.strand)
    }

    ## ... fracOverlap
    for (i in seq_along(fracOverlap)) {
        gr2 <- fracOverlap[[i]]
        if (is(gr2, "GRangesList")) {
            gr2 <- BiocGenerics::unlist(gr2)
        }
        ov <- GenomicRanges::findOverlaps(query = x, subject = gr2,
                                          ignore.strand = ignore.strand)
        tilesov <- pmin(GenomicRanges::end(x[S4Vectors::queryHits(ov)]),
                        GenomicRanges::end(gr2[S4Vectors::subjectHits(ov)])) -
            pmax(GenomicRanges::start(x[S4Vectors::queryHits(ov)]),
                 GenomicRanges::start(gr2[S4Vectors::subjectHits(ov)])) + 1
        tilesovSum <- tapply(tilesov, S4Vectors::queryHits(ov), sum)
        posPerTile <- rep(0, length(x))
        posPerTile[as.numeric(names(tilesovSum))] <- tilesovSum
        nm <- paste0(names(fracOverlap)[i], ".fracOverlap")
        df[[nm]] <- posPerTile / width(x)
    }

    ## ... numOverlap
    for (i in seq_along(numOverlap)) {
        gr2 <- numOverlap[[i]]
        nm1 <- paste0(names(numOverlap)[i], ".numOverlapWithin")
        nm2 <- paste0(names(numOverlap)[i], ".numOverlapAny")
        # remark: findOverlaps(..., type = "within") refers to the query
        ov1 <- GenomicRanges::findOverlaps(query = gr2, subject = x,
                                           type = "within",
                                           ignore.strand = ignore.strand)
        ov1perGr <- table(S4Vectors::subjectHits(ov1))
        nwithin <- rep(0L, length(x))
        nwithin[as.numeric(names(ov1perGr))] <- as.integer(ov1perGr)
        df[[nm1]] <- nwithin
        df[[nm2]] <- GenomicRanges::countOverlaps(query = x, subject = gr2,
                                                  type = "any",
                                                  ignore.strand = ignore.strand)
    }

    ## ... nearest
    for (i in seq_along(nearest)) {
        gr2 <- nearest[[i]]
        if (is(gr2, "GRangesList")) {
            gr2 <- BiocGenerics::unlist(gr2)
        }
        nm1 <- paste0(names(nearest)[i], ".nearestName")
        nm2 <- paste0(names(nearest)[i], ".nearestDistance")
        j <- GenomicRanges::nearest(x = x, subject = gr2, select = "arbitrary",
                                    ignore.strand = ignore.strand)
        df[[nm1]] <- names(gr2)[j]
        dists <- rep(NA_integer_, length(x))
        ok <- !is.na(j)
        dists[ok] <- GenomicRanges::distance(x = x[ok], y = gr2[j[ok]],
                                             ignore.strand = ignore.strand)
        df[[nm2]] <- dists
    }

    ## return result
    S4Vectors::mcols(x) <- df
    return(x)
}
