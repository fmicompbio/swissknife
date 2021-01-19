#' @title Calc mappability of a genome.
#'
#' @description Given a window width and the maximum number of allowed hits per
#'     window sequence, calculate the mappability of all overlapping windows in
#'     that genome.
#'
#' @author Michael Stadler
#'
#' @param genome The genome sequence to work on. Either a \code{\link[BSgenome]{BSgenome}}
#'     object, a \code{character} scalar with the name of an installed \code{\link[BSgenome]{BSgenome}}
#'     or with a file paht and name pointing to a fasta file with the genome sequence.
#' @param genomeIndex \code{character} scalar with the path to the bowtie index
#'     and prefix to align against, in the form \code{</path/to/index>/<prefix>},
#'     or the name of an installed \code{Rbowtie} index package created by the
#'     \pkg{QuasR} package for an installed \code{BSgenome} package.
#' @param kmerLength \code{numeric} scalar specifying the k-mer length (width of
#'     overlapping windows in \code{genome}), usually set to the typical read
#'     length for which to get the mappable regions.
#' @param maxHits \code{numeric} scalar specifying the maximum number of hits
#'     (matches) of a window sequence in \code{genome} to be considered mappable.
#' @param Ncpu \code{numeric} scalar specifying the number of CPU threads to use
#'     for alignments.
#' @param quiet \code{logical} scalar indicating if progress information should
#'     be printed on the console.
#'
#' @details Sequences of all overlapping windows are extracted from the genome
#'     and aligned back to that genome using \code{\link[Rbowtie]{bowtie}} with
#'     parameters \code{--concise -f -v 0 -a -B 1 -m maxHits}. If no more than
#'     \code{maxHits} hits are found, the window is defined mappable.
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object with mappable regions.
#'     All plus-strand sequences in \code{genome} of length \code{kmerLength}
#'     with their start (leftmost) position overlapping the \code{GRanges} object
#'     do not generate more than \code{maxHits} hits when aligned to the genome.
#'
#' @examples
#' library(Rbowtie)
#' 
#' genomefile <- system.file("extdata", "getMappableRegions", "hg19sub.fa", package = "swissknife")
#' indexdir <- tempfile()
#' indexpre <- "index"
#' indexname <- file.path(indexdir, indexpre)
#' idx <- bowtie_build(genomefile, indexdir)
#' 
#' mapgr <- getMappableRegions(genomefile, indexname, 50, quiet = FALSE)
#' mapgr
#'
#' @seealso \code{\link[Rbowtie]{bowtie}} in package \pkg{Rbowtie} used by
#'     \code{getMappableRegions} to align reads to the genome;
#'     \code{\link[Rbowtie]{bowtie_build}} in package \pkg{Rbowtie} for
#'     indexing a genome.
#'
#' @importFrom GenomicRanges GRanges gaps sort
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BSgenome seqinfo getSeq
#' @importFrom GenomeInfoDb seqnames seqlengths
#' @importFrom XVector width
#' 
#' @export
getMappableRegions <- function(genome,
                               genomeIndex,
                               kmerLength = 50,
                               maxHits = 1,
                               Ncpu = 2,
                               quiet = TRUE) {
    # check arguments
    # ... genome
    if(is(genome, "BSgenome")) {
        chrinfo <- BSgenome::seqinfo(genome)
        chrs <- BSgenome::getSeq(genome, names = seqnames(chrinfo),
                                 start = 1, end = seqlengths(chrinfo))
    } else if (is.character(genome) && length(genome) == 1L) {
        if (suppressWarnings(require(genome, character.only = TRUE, quietly = TRUE))) {
            genome <- get(genome)
            chrinfo <- BSgenome::seqinfo(genome)
            chrs <- BSgenome::getSeq(genome, names = seqnames(chrinfo),
                                     start = 1, end = seqlengths(chrinfo))
        } else if (file.exists(genome)) {
            chrs <- Biostrings::readDNAStringSet(genome)
        } else {
            stop("'genome' is neither a valid file nor a BSgenome object.")
        }
    }
    if (!quiet) {
        message("Will calculate mappability for ", length(chrs), " chromosomes ",
                "(", sum(width(chrs)), " base pairs)")
    }
    # ... genomeIndex
    if (is.character(genomeIndex) && length(genomeIndex) == 1L) {
        if (suppressWarnings(require(genomeIndex, character.only = TRUE, quietly = TRUE))) {
            genomeIndex <- file.path(system.file("alignmentIndex", package = genomeIndex), "bowtieIndex")
        } else {
            indexDir <- dirname(genomeIndex)
            indexPre <- basename(genomeIndex)
            if (!file.exists(indexDir)) {
                stop("index directory '", indexDir, "' does not exist.")
            } else if (!file.exists(file.path(indexDir, paste0(indexPre, ".1.ebwt")))) {
                stop("no bowtie index '", indexPre, "' found in '", indexDir, "'.")
            }
        }
    } else {
        stop("'genomeIndex' is not a valid bowtie index.")
    }
    # ... other arguments
    stopifnot(exprs = {
        is.numeric(kmerLength)
        length(kmerLength) == 1L
        kmerLength > 0
        is.numeric(maxHits)
        length(maxHits) == 1L
        maxHits > 0
        is.numeric(Ncpu)
        length(Ncpu) == 1L
        Ncpu > 0
    })
    
    # iterate over chromosomes
    nonmapL <- list()
    for (i in seq_along(chrs)) {
        if (!quiet)
            message("    ", names(chrs)[i], "...", appendLF = FALSE)

        # write overlapping window sequences to a file
        readfile <- .writeWindowsToTempFile(as.character(chrs[[i]]), kmerLength)

        # align to the genome
        outfiles <- .alignWindowsToGenome(readfile, genomeIndex, maxHits, Ncpu)

        # identify mapped / unmapped windows
        idL <- lapply(outfiles[c("fun", "fmax")], function(f) {
            ids <- numeric(0)
            if (file.exists(f)) {
                tmp <- readLines(f)
                stopifnot(length(tmp) %% 2 == 0)
                ids <- as.numeric(sub("^>", "", tmp[seq(1, length(tmp) - 1, by = 2)]))
            }
            ids
        })
        nonmapL[[i]] <- GRanges(names(chrs)[i], IRanges(start = unlist(idL), width = 1),
                                seqlengths = seqlengths(chrs), strand = "+")
        ## remark: this will also lead to the last kmerLength-1 bases as mappable

        # clean up
        unlink(readfile)
        unlink(outfiles)
        if (!quiet)
            message("done (",
                    round(length(nonmapL[[i]]) * 100 / (width(chrs)[i] - kmerLength + 1), 1),
                    "% non-mappable)")
    }
    
    # return result
    gr <- gaps(do.call(c, nonmapL))
    gr <- sort(gr[strand(gr) == "+"])
    if (!quiet)
        message("Total mappability: ",
                round(sum(width(gr)) * 100 / sum(width(chrs) - kmerLength + 1), 1),
                "% in ", length(gr), " segments")
    return(gr)
}

# write all sub-strings of length 'w' in 'chr' to 'fname'
.writeWindowsToTempFile <- function(chr, w, fname = NULL) {
    if (is.null(fname))
        fname <- tempfile(fileext = ".fa")
    s <- seq.int(nchar(chr) - w + 1)
    wseqs <- paste0(">", as.character(s), "\n",
                    substr(rep(chr, length(s)), start = s, stop = s + w - 1))
    fh <- file(fname, open = "wt")
    writeLines(wseqs, con = fh)
    close(fh)
    return(fname)
}

# align 'fname' to 'index'
.alignWindowsToGenome <- function(fname, index, m = 1, p = 1,
                                  fmax = NULL, fun = NULL, fout = NULL) {
    if (is.null(fmax))
        fmax <- tempfile()
    if (is.null(fun))
        fun <- tempfile()
    if (is.null(fout))
        fout <- tempfile()
    args <- sprintf("--max %s --un %s -f -p %d -v 0 -a -B 1 --quiet -m %d %s %s %s",
                    fmax, fun, p, m, index, fname, fout)
    res <- system2(command = system.file("bowtie", package = "Rbowtie"), args = args,
                   stdout = FALSE, stderr = FALSE)
    return(c(fmax = fmax, fun = fun, fout = fout))
}