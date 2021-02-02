#' @title Get regions tiling a genome.
#'
#' @description Get sequential, potentially annotated regions of a fixed lengths
#'    (tiles) along chromosomes of a genome.
#'
#' @author Michael Stadler
#'
#' @param genome The genome to work on. Either a \code{\link[BSgenome]{BSgenome}}
#'     object, a \code{character} scalar with the name of an installed \code{\link[BSgenome]{BSgenome}}
#'     or with a file path and name pointing to a fasta file with the genome sequence,
#'     or a named \code{numeric} vector giving the names and lengths of chromosomes.
#' @param tileWidth \code{numeric} scalar with the tile length.
#' @param hasOverlap Named \code{list} with \code{\link[GenomicRanges]{GRanges}} object(s).
#'     For each list element, a logical vector "X.hasOverlap" will be added to the
#'     \code{mcols} of the result, with \code{TRUE} for each tile that overlaps
#'     any region in that element. "X" is obtained from \code{names(hasOverlap)}.
#' @param fracOverlap  Named \code{list} with \code{\link[GenomicRanges]{GRanges}} object(s).
#'     For each list element, a numeric vector "X.fracOverlap" will be added to the
#'     \code{mcols} of the result, with a value between 0 and 1 giving the fraction
#'     of bases in a tile that overlaps with any region in that element. "X" is
#'     obtained from \code{names(fracOverlap)}.
#' @param numOverlap  Named \code{list} with \code{\link[GenomicRanges]{GRanges}} object(s).
#'     For each list element, two numeric vectors "X.numOverlapWithin" and
#'     "X.numOverlapAny" will be added to the \code{mcols} of the result, giving
#'     the number of ranges in that element that are fully contained within
#'     a tile, or that overlap with a tile in any way, respectively. "X" is
#'     obtained from \code{names(numOverlap)}.
#' @param nearest Named \code{list} with \code{\link[GenomicRanges]{GRanges}} object(s).
#'     For each list element, two numeric vectors "X.nearestName" and
#'     "X.nearestDistance" will be added to the \code{mcols} of the result, giving
#'     the name and distance of the nearest range in that element for each tile. "X" is
#'     obtained from \code{names(nearest)}, and the values of "X.nearestName" from
#'     \code{names(nearest$X)}.
#'
#' @details The last tile in each chromosome is dropped if it would be shorter
#'     than \code{tileWidth}. Generated tiles are unstranded (\code{*}) and
#'     therefore overlaps or searching for nearest neighbors are ignoring
#'     strands of annotations (\code{ignore.strand=TRUE}). If multiple nearest
#'     ranges are at the same distance from a tile, an arbitrary one is
#'     reported in "X.nearestName".
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object with genome tiling regions.
#'     Optional tile annotations are contained in its metadata columns (\code{mcols}).
#'
#' @examples
#' library(GenomicRanges)
#' 
#' tss <- GRanges("chr1", IRanges(c(1, 10, 30), width = 1,
#'                                names = paste0("t", 1:3)))
#' blacklist <- GRanges("chr1", IRanges(20, width = 5))
#' getGenomicTiles(c(chr1 = 45, chr2 = 12), tileWidth = 10,
#'                 hasOverlap = list(Blacklist = blacklist),
#'                 fracOverlap = list(Blacklist = blacklist),
#'                 numOverlap = list(TSS = tss),
#'                 nearest = list(TSS = tss))
#'
#' @seealso \code{\link[GenomicRanges]{tileGenome}}, \code{\link[GenomicRanges]{findOverlaps}}
#'     and \code{\link[GenomicRanges]{nearest}} in package \pkg{GenomicRanges} used by
#'     \code{getGenomicTiles} internally.
#'
#' @importFrom GenomicRanges GRanges tileGenome width findOverlaps nearest distance
#' @importFrom BSgenome seqinfo
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Biostrings fasta.seqlengths
#' @importFrom methods is
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @export
getGenomicTiles <- function(genome,
                            tileWidth,
                            hasOverlap = list(),
                            fracOverlap = list(),
                            numOverlap = list(),
                            nearest = list()) {
    ## check arguments
    ## ... genome
    if(is(genome, "BSgenome")) {
        chrlens <- GenomeInfoDb::seqlengths(BSgenome::seqinfo(genome))
    } else if (is.character(genome) && length(genome) == 1L) {
        if (suppressWarnings(require(genome, character.only = TRUE, quietly = TRUE))) {
            genome <- get(genome)
            chrlens <- GenomeInfoDb::seqlengths(BSgenome::seqinfo(genome))
        } else if (file.exists(genome)) {
            chrlens <- Biostrings::fasta.seqlengths(genome)
        } else {
            stop("'genome' is neither a valid file nor a BSgenome object.")
        }
    } else if (is.numeric(genome) && !is.null(names(genome))) {
        chrlens <- genome
    } else {
        stop("'genome' is not a valid argument for getGenomicTiles()")
    }
    ## ... other arguments
    stopifnot(exprs = {
        # tileWidth
        is.numeric(tileWidth)
        length(tileWidth) == 1L
        tileWidth > 0
        # hasOverlap
        is.list(hasOverlap)
        length(hasOverlap) == 0L || !is.null(names(hasOverlap))
        all(unlist(lapply(hasOverlap, is, "GRanges")))
        # fracOverlap
        is.list(fracOverlap)
        length(fracOverlap) == 0L || !is.null(names(fracOverlap))
        all(unlist(lapply(fracOverlap, is, "GRanges")))
        # numOverlap
        is.list(numOverlap)
        length(numOverlap) == 0L || !is.null(names(numOverlap))
        all(unlist(lapply(numOverlap, is, "GRanges")))
        # nearest
        is.list(nearest)
        length(nearest) == 0L || !is.null(names(nearest))
        all(unlist(lapply(nearest, is, "GRanges")))
        all(unlist(lapply(nearest, function(x) !is.null(names(x)))))
    })
    
    ## create tiles
    gr <- GenomicRanges::tileGenome(seqlengths = chrlens, tilewidth = tileWidth,
                                    cut.last.tile.in.chrom = TRUE)
    gr <- gr[GenomicRanges::width(gr) == tileWidth]
    
    ## annotate tiles
    df <- S4Vectors::mcols(gr)
    
    ## ... hasOverlap
    for (i in seq_along(hasOverlap)) {
        nm <- paste0(names(hasOverlap)[i], ".hasOverlap")
        df[[nm]] <- IRanges::overlapsAny(query = gr, subject = hasOverlap[[i]],
                                         ignore.strand = TRUE)
    }

    ## ... fracOverlap
    for (i in seq_along(fracOverlap)) {
        gr2 <- fracOverlap[[i]]
        ov <- GenomicRanges::findOverlaps(query = gr, subject = gr2,
                                          ignore.strand = TRUE)
        tilesov <- pmin(GenomicRanges::end(gr[S4Vectors::queryHits(ov)]),
                        GenomicRanges::end(gr2[S4Vectors::subjectHits(ov)])) -
            pmax(GenomicRanges::start(gr[S4Vectors::queryHits(ov)]),
                 GenomicRanges::start(gr2[S4Vectors::subjectHits(ov)])) + 1
        tilesovSum <- tapply(tilesov, S4Vectors::queryHits(ov), sum)
        posPerTile <- rep(0, length(gr))
        posPerTile[as.numeric(names(tilesovSum))] <- tilesovSum
        nm <- paste0(names(fracOverlap)[i], ".fracOverlap")
        df[[nm]] <- posPerTile / tileWidth
    }

    ## ... numOverlap
    for (i in seq_along(numOverlap)) {
        gr2 <- numOverlap[[i]]
        nm1 <- paste0(names(numOverlap)[i], ".numOverlapWithin")
        nm2 <- paste0(names(numOverlap)[i], ".numOverlapAny")
        # remark: findOverlaps(..., type = "within") refers to the query
        ov1 <- GenomicRanges::findOverlaps(query = gr2, subject = gr,
                                           type = "within",
                                           ignore.strand = TRUE)
        ov1perGr <- table(S4Vectors::subjectHits(ov1))
        nwithin <- rep(0L, length(gr))
        nwithin[as.numeric(names(ov1perGr))] <- as.integer(ov1perGr)
        df[[nm1]] <- nwithin
        df[[nm2]] <- GenomicRanges::countOverlaps(query = gr, subject = gr2,
                                                  type = "any",
                                                  ignore.strand = TRUE)
    }

    ## ... nearest
    for (i in seq_along(nearest)) {
        gr2 <- nearest[[i]]
        nm1 <- paste0(names(nearest)[i], ".nearestName")
        nm2 <- paste0(names(nearest)[i], ".nearestDistance")
        j <- GenomicRanges::nearest(x = gr, subject = gr2, select = "arbitrary",
                                    ignore.strand = TRUE)
        df[[nm1]] <- names(gr2)[j]
        dists <- rep(NA_integer_, length(gr))
        ok <- !is.na(j)
        dists[ok] <- GenomicRanges::distance(x = gr[ok], y = gr2[j[ok]],
                                             ignore.strand = TRUE)
        df[[nm2]] <- dists
    }

    ## return result
    S4Vectors::mcols(gr) <- df
    return(gr)
}
