#' @title Tabulate insert sizes from paired-end alignments in bam files.
#'
#' @description Read and tabulate the insert sizes from paired-end alignments
#'   contained in one or several bam files. By default, all properly aligned
#'   read pairs are included. Optionally, alignments can be restricted to
#'   those in a specific genomic region region (\code{regions} argument) or
#'   the number of alignments read can be limited (\code{nmax} argument).
#'
#' @author Michael Stadler
#'
#' @param fname \code{character} vector with one or several bam files. If
#'   multiple files are given, insert sizes from all will be pooled and
#'   tabulated together.
#' @param regions \code{GRanges} object. Only alignments falling into these
#'   regions will be used. If \code{NULL} (the default), all alignments are
#'   used.
#' @param nmax \code{numeric(1)} specifying the maximal number of alignments to
#'   read. If \code{NA} (the default), the alignments in \code{regions} (if
#'   \code{regions} are not \code{NULL}) or in the bam file will be used.
#' @param isizemax \code{numeric(1)} specifying the maximal insert size to
#'   report. Larger insert sizes will be set to \code{isizemax} with a warning.
#'
#' @return \code{integer} vector with the number of insert sizes. The element at
#'   position \code{i} gives the observed number of alignment pairs with an
#'   insert size of \code{i}.
#'
#' @seealso \code{\link[Rsamtools]{scanBam}} used to read alignments.
#'
#' @examples
#' if (requireNamespace("Rsamtools", quietly = TRUE)) {
#'     bamf <- system.file("extdata", "getInsertSizeDistFromBam", "atac_mm10.bam",
#'                         package = "swissknife")
#'     isize <- getInsertSizeDistFromBam(bamf)
#'     plot(isize, type = "l",
#'          xlab = "Insert size (bp)", ylab = "Number of fragments")
#' }
#'
#' @importFrom IRanges IRangesList
#' 
#' @export
getInsertSizeDistFromBam <- function(fname,
                                     regions = NULL,
                                     nmax = NA_integer_,
                                     isizemax = 800) {
    .assertPackagesAvailable(c("Rsamtools"))
    
    # set scanning parameters
    if (is.null(regions)) {
        regions <- IRanges::IRangesList()
    }
    params <- Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(isPaired = TRUE,
                                      isProperPair = TRUE,
                                      isUnmappedQuery = FALSE,
                                      hasUnmappedMate = FALSE,
                                      isFirstMateRead = TRUE,
                                      isSecondaryAlignment = FALSE,
                                      isNotPassingQualityControls = NA,
                                      isDuplicate = FALSE,
                                      isSupplementaryAlignment = FALSE),
        what = c("isize"),
        which = regions)
    
    # open bam file(s)
    bfl <- Rsamtools::BamFileList(fname, yieldSize = min(2e6, nmax, na.rm = TRUE))
    
    # read insert sizes
    bfl <- open(bfl)
    isizeL <- lapply(bfl, function(x) {
        itab <- integer(isizemax)
        while (length(i1 <- unlist(Rsamtools::scanBam(x, param = params), use.names = FALSE)) > 0) {
            itab <- itab + tabulate(abs(i1), isizemax)
        }
        itab
    })
    bfl <- close(bfl)
    
    
    # combine and return results
    isize <- Reduce(f = "+", x = isizeL, init = integer(isizemax))
    return(isize)
}
