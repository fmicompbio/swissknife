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
#' @param outFile Character scalar, giving the path to which the script 
#'     should be copied. The path is relative to the root of the active 
#'     project. If a file with this name already exists, the function will 
#'     ask for confirmation before overwriting it. 
#' @param copyTests Logical scalar, defining whether to copy unit tests for 
#'     the utility functions to \code{tests/testthat/tests-<name>.R}, 
#'     where \code{<name>} is the base name of \code{outFile}.
#'     
#' @author Charlotte Soneson
#' 
#' @export
#' 
#' @importFrom usethis use_template use_testthat
addUtilsFunctions <- function(outFile = "R/utils.R", copyTests = TRUE) {
    usethis::use_template("utils_template.R", package = "swissknife",
                          save_as = outFile)
    if (copyTests) {
        usethis::use_testthat()
        usethis::use_template("tests_utils_template.R", package = "swissknife",
                              save_as = file.path("tests", "testthat", 
                                                  paste0("tests-", 
                                                         basename(outFile))))
    }
}
