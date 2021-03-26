#' Make sure that packages are available
#' 
#' Try loading the namespaces of the packages given in \code{pkgs}, and
#' throw an exception with an informative error message if that is not
#' the case.
#' 
#' @param pkgs Character vector with package names.
#' @param bioc Logical scalar. If \code{TRUE}, suggest to install the missing
#'   package(s) using \code{BiocManager} instead of \code{install.packages}.
#' 
#' @author Michael Stadler
#' 
.assertPackagesAvailable <- function(pkgs, bioc = TRUE) {
    stopifnot(exprs = {
        is.character(pkgs)
        is.logical(bioc)
        length(bioc) == 1L
    })
    
    avail <- unlist(lapply(pkgs,
                           function(pkg) {
                               requireNamespace(pkg, quietly = TRUE)
                           }))
    
    if (any(!avail)) {
        caller <- deparse(sys.calls()[[sys.nframe() - 1]])
        callerfunc <- sub("\\(.+$", "", caller)
        stop("The package", ifelse(sum(!avail) > 1, "s '", " '"),
             paste(pkgs[!avail], collapse = "', '"), "' ",
             ifelse(sum(!avail) > 1, "are", "is"), " required for ",
             callerfunc, "(), but not installed.\n",
             "Install ", ifelse(sum(!avail) > 1, "them", "it"), " using:\n",
             ifelse(bioc, "BiocManager::install(c(\"", "install.packages(c(\""),
             paste(pkgs[!avail], collapse = "\", \""), "\"))", call. = FALSE)
    }
    
    invisible(TRUE)
}
