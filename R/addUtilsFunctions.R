#' Copy utility functions to package
#' 
#' This function copies handy utility functions to a new script in a 
#' specified location. Currently, the script contains the following utility 
#' functions:
#' \itemize{
#' \item \code{.assertScalar()} - convenience function to check the validity 
#'     of scalar variables.
#' \item \code{.assertVector()} - convenience function to check the validity 
#'     of vector variables. 
#' }
#' 
#' @param outfile Character scalar, giving the path to which the script 
#'     should be copied. The path is relative to the root of the active 
#'     project. If a file with this name already exists, the function will 
#'     ask for confirmation before overwriting it. 
#'     
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @importFrom usethis use_template
addUtilsFunctions <- function(outfile = "R/utils.R") {
    usethis::use_template("utils_template.R", package = "swissknife",
                          save_as = outfile)
}