#' @title Parse R and R package versions from session informations
#'
#' @author Michael Stadler
#'
#' @description The function parses the R version and R package versions from
#'   session information (created by \code{sessionInfo()}, tested with R 3.6) in
#'   files provided in \code{infiles}. Two types of files are currently
#'   supported: \itemize{ \item{Rout: }{Files containing R console output
#'   (created by \code{R CMD BATCH scrip.R output.Rout}} \item{md: }{Files
#'   containing markdown output created by \code{rmarkdown::render('input.Rmd',
#'   clean = FALSE)}, which will keep the intermediate \code{.md} file.}}
#'
#' @param infiles Character vector with text files (extension must be either
#'   \code{.Rout} or \code{.md}), containing session information to parse out.
#'
#' @return A \code{list} of \code{list}s with one element in the outer list for
#'   each R version, contianing an innter list with elements \code{files} and
#'   \code{packages}.
#'
#' @examples
#' f <- list.files(system.file("extdata", "parsePkgVersions",
#'                             package = "swissknife"),
#'                 full.names = TRUE)
#' parsePkgVersions(f)
#'
#' @importFrom utils sessionInfo
#' @importFrom tools file_ext
#'
#' @export
#' 
parsePkgVersions <- function(infiles) {
    
    ## pre-flight checks
    stopifnot(exprs = {
        is.character(infiles)
        all(file.exists(infiles))
        all(tools::file_ext(infiles) %in% c("Rout", "md"))
    })
    
    ## parse files
    all_packages <- list()
    fext <- tools::file_ext(infiles)
    pat1 <- c(Rout = "> sessionInfo()",
              md = "sessionInfo()")[fext]
    pat2 <- c(Rout = "other attached packages:",
              md = "## other attached packages:")[fext]
    pat3 <- c(Rout = "loaded via a namespace (and not attached):",
              md = "## loaded via a namespace (and not attached):")[fext]
    pat4 <- c(Rout = "\\[",
              md = "^\\[[0-9]+\\]$")[fext]
    getPkgFun <- c(Rout = function(s) setdiff(setdiff(strsplit(s, " ")[[1]], " "), ""),
                   md = function(s) setdiff(strsplit(s, " ")[[1]], c("", " ", "##")))[fext]

    for (i in seq_along(infiles)) {
        x <- readLines(infiles[i])
        idx1 <- which(x == pat1[i])
        idx2 <- which(x == pat2[i])
        idx3 <- which(x == pat3[i])
        
        if (length(idx1) == 1L & length(idx2) == 1L & length(idx3) == 1L) {
            thisR <- switch(fext[i], Rout = x[idx1 + 1], md = sub("^## ","",x[idx1 + 4]))

            if (!(thisR %in% names(all_packages))) {
                all_packages[[thisR]] <- list(files = character(0), packages = character(0))
            }
            
            all_packages[[thisR]][["files"]] <- c(all_packages[[thisR]]$files, infiles[i])
            all_packages[[thisR]][["packages"]] <-
                unique(c(all_packages[[thisR]][["packages"]],
                         do.call(c, lapply((idx2 + 1):(idx3 - 2), function(j) {
                             grep(pat4[i],
                                  getPkgFun[[i]](x[j]),
                                  value = TRUE,
                                  invert = TRUE)
                         }))))
        } else {
            warning("skipping ", infiles[i], " (no or multiple session informations found)")
        }
    }
    
    ## sort
    all_packages <- lapply(all_packages, function(x) {
        x$files <- sort(x$files)
        x$packages <- sort(x$packages)
        x
    })
    if (length(all_packages) > 0) {
        all_packages <- all_packages[ order(names(all_packages)) ]
    }
    
    ## return
    all_packages
}

