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
#' @param hasOverlap,fracOverlap,numOverlap,nearest Named \code{list}s with
#'     \code{\link[GenomicRanges]{GRanges}} or \code{\link[GenomicRanges]{GRangesList}}
#'     object(s) used to annotate genomic tiles. See \code{\link{annotateRegions}}
#'     for details.
#' @param addSeqComp \code{logical} scalar. If \code{TRUE} and primary sequence
#'     can be obtained from \code{genome}, also add sequence composition features
#'     for each tile to the annotations. Currently, the following features are
#'     included: percent of G+C bases ("percGC"), CpG observed-over-expected ratio
#'     ("CpGoe").
#'
#' @details The last tile in each chromosome is dropped if it would be shorter
#'     than \code{tileWidth}. Generated tiles are unstranded (\code{*}) and
#'     therefore overlaps or searching for nearest neighbors are ignoring
#'     strands of annotations (\code{ignore.strand=TRUE}).
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
#' @seealso \code{\link[GenomicRanges]{tileGenome}} and \code{\link{annotateRegions}}
#'     used by \code{getGenomicTiles} internally.
#'
#' @importFrom GenomicRanges GRanges tileGenome width findOverlaps nearest distance
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Biostrings fasta.seqlengths readDNAStringSet oligonucleotideFrequency
#' @importFrom methods is
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics unlist
#' 
#' @export
getGenomicTiles <- function(genome,
                            tileWidth,
                            hasOverlap = list(),
                            fracOverlap = list(),
                            numOverlap = list(),
                            nearest = list(),
                            addSeqComp = TRUE) {
    ## check arguments
    .assertPackagesAvailable("BSgenome")
    ## ... genome and addSeqComp
    stopifnot(exprs = {
        is.logical(addSeqComp)
        length(addSeqComp) == 1L
    })
    if (is(genome, "BSgenome")) {
        chrlens <- GenomeInfoDb::seqlengths(BSgenome::seqinfo(genome))
        genomeobj <- genome
    } else if (is.character(genome) && length(genome) == 1L) {
        if (suppressWarnings(require(genome, character.only = TRUE, quietly = TRUE))) {
            genomeobj <- get(genome)
            chrlens <- GenomeInfoDb::seqlengths(BSgenome::seqinfo(genomeobj))
        } else if (file.exists(genome)) {
            chrlens <- Biostrings::fasta.seqlengths(genome)
            genomeobj <- Biostrings::readDNAStringSet(genome)
        } else {
            stop("'genome' is neither a valid file nor a BSgenome object.")
        }
    } else if (is.numeric(genome) && !is.null(names(genome))) {
        chrlens <- genome
        if (addSeqComp) {
            warning("ignoring 'addSeqComp' (sequence not available from 'genome')")
            addSeqComp <- FALSE
        }
    } else {
        stop("'genome' is not a valid argument for getGenomicTiles()")
    }
    ## ... other arguments
    stopifnot(exprs = {
        # tileWidth
        is.numeric(tileWidth)
        length(tileWidth) == 1L
        tileWidth > 0
    })
    
    ## create tiles
    gr <- GenomicRanges::tileGenome(seqlengths = chrlens, tilewidth = tileWidth,
                                    cut.last.tile.in.chrom = TRUE)
    gr <- gr[GenomicRanges::width(gr) == tileWidth]
    
    ## annotate tiles
    ## ... addSeqComp
    if (addSeqComp) {
        df <- S4Vectors::mcols(gr)
        grSeq <- BSgenome::getSeq(x = genomeobj, names = gr)
        grF1 <- Biostrings::oligonucleotideFrequency(x = grSeq, width = 1L,
                                                     as.prob = TRUE,
                                                     simplify.as = "matrix")
        grF2 <- Biostrings::oligonucleotideFrequency(x = grSeq, width = 2L,
                                                     as.prob = TRUE,
                                                     simplify.as = "matrix")
        df[["percGC"]] <- rowSums(grF1[,c("C","G")]) * 100
        df[["CpGoe"]] <- grF2[, "CG"] / (grF1[, "C"] * grF1[, "G"])
        S4Vectors::mcols(gr) <- df
    }
    ## ... overlap and nearest annotations
    gr <- annotateRegions(x = gr, hasOverlap = hasOverlap,
                          fracOverlap = fracOverlap, numOverlap = numOverlap,
                          nearest = nearest, ignore.strand = TRUE)

    ## return result
    return(gr)
}
