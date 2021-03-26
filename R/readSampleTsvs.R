utils::globalVariables(c("V1","V2","SAMPLE_ID")) # prevent R CMD check from complaining

#' @title Read sample tsv files from seqdata storage
#' 
#' @author Charlotte Soneson
#' 
#' @description The function searches the provided \code{seqdataDir} for tsv
#'   files corresponding to the provided \code{sampleIds} and returns a
#'   \code{data.frame} containing the metadata for all these samples.
#'   
#' @param seqdataDir Character scalar, the path to the directory containing the
#'   tsv files.
#' @param sampleIds Character vector with sample IDs, which will be matched
#'   against the file names in \code{seqDataDir}. The sample IDs should not
#'   contain the \code{.tsv} suffix.
#' @param keepMulti Logical scalar, indicating whether to keep samples that
#'   match more than one tsv file. If \code{TRUE}, these samples are represented
#'   by multiple rows in the table. If \code{FALSE}, these samples are excluded.
#'   In any case, a warning will be generated, listing the samples with multiple
#'   matching files.
#' @param ... Additional arguments that will be passed to \code{list.files},
#'   e.g. to make the search case-insensitive or search recursively.
#' 
#' @return A \code{data.frame} with metadata for the provided \code{sampleIds}.
#' 
#' @examples
#' if (requireNamespace("dplyr") && requireNamespace("tidyr")) {
#'     print(readSampleTsvs(seqdataDir = system.file("extdata/readSampleTsvs", 
#'                                                   package = "swissknife"), 
#'                          sampleIds = c("readSampleTsvsEx1",
#'                                        "readSampleTsvsEx2",
#'                                        "readSampleTsvsEx3")))
#' }
#' 
#' @importFrom utils read.delim
#' 
#' @export
#' 
readSampleTsvs <- function(seqdataDir = "/tungstenfs/groups/gbioinfo/seqdata", 
                           sampleIds, keepMulti = TRUE, ...) {
    
    ## Check if dplyr and tidyr are available
    .assertPackagesAvailable(c("dplyr", "tidyr"), bioc = FALSE)
    
    ## List all tsv files in seqdataDir matching any of the sample IDs
    matchingFiles <- list.files(path = seqdataDir, 
                                pattern = paste(paste0(sampleIds, ".*\\.tsv"),
                                                collapse = "|"),
                                full.names = TRUE, ...)
    
    ## Go through provided sample IDs
    do.call(dplyr::bind_rows, lapply(sampleIds, function(s) {
        ## Get the matching file(s) for sample s
        f <- grep(paste0(s, ".*\\.tsv"), matchingFiles,
                  value = TRUE)

        if (length(f) == 0) {
            ## If no matching file is found, the sample is excluded
            warning("No file matching sample ID ", s, " in ", 
                    seqdataDir, 
                    ".\nThis sample ID is not included in the final table.",
                    call. = FALSE)
            data.frame()
        } else if (length(f) > 1 && !keepMulti) {
            ## If multiple matching files are found and keepMulti is FALSE, the
            ## sample is excluded.
            warning("More than one file matching sample ID ", s, " in ", 
                    seqdataDir, ":\n",
                    paste(basename(f), collapse = "\n"),
                    ".\nThis sample ID is not included in the final table.", 
                    call. = FALSE)
            data.frame()
        } else {
            if (length(f) > 1) {
                ## If multiple matching files are found and keepMulti is TRUE,
                ## generate a warning, but continue.
                warning("More than one file matching sample ID ", s, " in ", 
                        seqdataDir, ":\n",
                        paste(basename(f), collapse = "\n"),
                        ".\nThis sample ID is represented by multiple lines ", 
                        "in the final table.", 
                        call. = FALSE)
            }
            ## Read matching file(s) and generate data.frame
            do.call(dplyr::bind_rows, lapply(f, function(f0) {
                utils::read.delim(f0, header = FALSE, as.is = TRUE) %>%
                    tidyr::spread(V1, V2) %>%
                    dplyr::mutate(SAMPLE_ID = s,
                                  TSV_FILE = basename(f0)) %>%
                    dplyr::select(SAMPLE_ID, dplyr::everything())
            }))
        }
    }))
}
