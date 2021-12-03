#' Utility function that makes sure that packages are available
#' 
#' The function tries loading the namespaces of the packages given in
#' \code{pkgs}, and throws an exception with an informative error message if
#' that is not the case.
#' 
#' @param pkgs Character vector with package names. Can be either just a
#'   package name or a string of the form \code{"githubuser/packagename"} for
#'   packages hosted on GitHub.
#' @param suggestInstallation Logical scalar. If \code{TRUE}, include an
#'   expression to install the missing package(s) as part of the generated
#'   error message.
#' 
#' @author Michael Stadler, Charlotte Soneson
#' 
#' @noRd
#' @keywords internal
.assertPackagesAvailable <- function(pkgs, suggestInstallation = TRUE) {
    stopifnot(exprs = {
        is.character(pkgs)
        is.logical(suggestInstallation)
        length(suggestInstallation) == 1L
    })
    
    avail <- unlist(lapply(sub("^[^/]+/", "", pkgs),
                           function(pkg) {
                               requireNamespace(pkg, quietly = TRUE)
                           }))
    
    if (any(!avail)) {
        caller <- deparse(sys.calls()[[sys.nframe() - 1]])
        callerfunc <- sub("\\(.+$", "", caller)
        haveBioc <- requireNamespace("BiocManager", quietly = TRUE)
        msg <- paste0("The package", ifelse(sum(!avail) > 1, "s '", " '"),
                      paste(sub("^[^/]+/", "", pkgs[!avail]), collapse = "', '"),
                      "' ",
                      ifelse(sum(!avail) > 1, "are", "is"), " required for ",
                      callerfunc, "(), but not installed.\n")
        if (suggestInstallation) {
            msg <- paste0(msg,
                          "Install ", ifelse(sum(!avail) > 1, "them", "it"), " using:\n",
                          ifelse(haveBioc, "", "install.packages(\"BiocManager\")\n"),
                          "BiocManager::install(c(\"",
                          paste(pkgs[!avail], collapse = "\", \""), "\"))")
        }
        stop(msg, call. = FALSE)
    }
    
    invisible(TRUE)
}
