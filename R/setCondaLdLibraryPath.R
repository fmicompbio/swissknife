#' Set LD_LIBRARY_PATH in conda environment
#' 
#' Create two shell scripts in subfolders of a given conda environment, 
#' setting and unsetting the LD_LIBRARY_PATH environment variable when 
#' activating and deactivating the environment. 
#' 
#' See https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables
#' for more information about setting environment variables in conda 
#' environments. To get the environment path for a \code{basilisk} package, 
#' use e.g. \code{basilisk::obtainEnvironmentPath()}. 
#' 
#' @param envPath Path to a conda environment for which the 
#'     \code{LD_LIBRARY_PATH} variable will be set. 
#' @param condaLibPath Path to the conda library folder that will be prepended 
#'     to the \code{LD_LIBRARY_PATH} environment variable. If \code{NULL} (the 
#'     default), it will be set to the subdirectory \code{lib/} of the 
#'     \code{basilisk}-managed conda instance, which is suitable for setting 
#'     the \code{LD_LIBRARY_PATH} for \code{basilisk} packages.
#'     For other applications, please set the conda path accordingly. The 
#'     path should point to a folder containing the C++ standard libraries. 
#' 
#' @author Charlotte Soneson
#' 
#' @export
#' @return Nothing is returned. The function is called for its side effect 
#'     of creating shell scripts in the provided conda environment.
#'     
#' @importFrom basilisk.utils getCondaDir
#' 
setCondaLdLibraryPath <- function(envPath, condaLibPath = NULL) {
    if (length(envPath) != 1 || !is.character(envPath)) {
        stop("envPath must be a character scalar")
    }
    if (!is.null(condaLibPath) && 
        (length(condaLibPath) != 1 || !is.character(condaLibPath))) {
        stop("condaLibPath must be a character scalar or NULL")
    }
    
    if (is.null(condaLibPath)) {
        condaLibPath <- file.path(basilisk.utils::getCondaDir(), "lib")
    }
    if (!file.exists(condaLibPath)) {
        stop("The provided conda path doesn't exist: ", condaLibPath)
    }
    ## Create subdirectories of the conda environment
    dir.create(file.path(envPath, "etc", "conda", "activate.d"), 
               showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(envPath, "etc", "conda", "deactivate.d"), 
               showWarnings = FALSE, recursive = TRUE)

    ## Write files
    writeLines(
        paste0("#!/bin/sh", 
               "\n\n", 
               "export TEMP_OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH", 
               "\n",
               "export LD_LIBRARY_PATH=", condaLibPath, ":$LD_LIBRARY_PATH"), 
        con = file.path(envPath, "etc", "conda", "activate.d", "env_vars.sh")
    )

    writeLines(
        paste0("#!/bin/sh", 
               "\n\n", 
               "export LD_LIBRARY_PATH=$TEMP_OLD_LD_LIBRARY_PATH",
               "\n",
               "unset TEMP_OLD_LD_LIBRARY_PATH"), 
        con = file.path(envPath, "etc", "conda", "deactivate.d", "env_vars.sh")
    )
}