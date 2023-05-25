#' @title Access example data
#'
#' @description Make example data available, typically for use in teaching.
#'
#' @author Michael Stadler
#'
#' @param name An optional character scalar specifying the data set(s) to
#'     be made available. The special name \code{"list"} (default) is used to
#'     print a data frame of available data sets with descriptions. The special
#'     name \code{"new"} will select the newest data set(s) available.
#' @param envir specifies the environment in which the data should be made
#'     available. By default, \code{envir = globalenv()}, which creates the
#'     example data objects in the user workspace. Possible alternative
#'     environment are for example \code{parent.frame()}, which is the
#'     environment in which \code{exampleData()} was called.
#' @param verbose A logical scalar. If \code{TRUE}, report what is being
#'     selected and made available.
#'
#' @return A \code{data.frame} (invisibly) with one row for each dataset
#'     that was made available in the global environment.
#'
#' @examples
#' exampleData()
#' exampleData("mycars")
#'
#' @import datasets
#'
#' @export
exampleData <- function(name = "list", envir = globalenv(), verbose = TRUE) {
    # internal data.frame with available packages
    .info <- data.frame(name = c("mycars",
                                 "mycars_summary",
                                 "2023-05"),
                        created = c("2023-05-23",
                                    "2023-05-23",
                                    "2023-05-24"),
                        description = c("re-encoded version of `datasets::mtcars`",
                                        "summary of `datasets::mtcars` (`hp` by `cyl`)",
                                        "collection: mycars, mycars_summary"),
                        is_collection = c(FALSE,
                                          FALSE,
                                          TRUE))

    # Check if dplyr and rlang are available
    .assertPackagesAvailable(c("dplyr", "rlang"))

    # check arguments
    if (!is.character(name) || length(name) != 1L) {
        stop("`name` must be a character scalar")
    }
    if (! name %in% c("list", "new", .info$name)) {
        stop("`name` must be one of: ", paste(.info$name, collapse = ", "))
    }
    if (! is.environment(envir)) {
        stop("`envir` must be an environment")
    }
    if (!is.logical(verbose) || length(verbose) != 1L) {
        stop("`verbose` must be a logical scalar")
    }

    # list or load data
    df <- switch (name,
        "list" = {
            print(.info)
            .info
        },
        "new" = {
            i <- which.max(as.Date(.info$created))
            if (verbose) {
                message("selecting newest dataset(s) `", .info[i, "name"], "`")
            }
            df <- exampleData(.info[i, "name"], envir, verbose)
            df
        },
        "mycars" = {
            mycars <- datasets::mtcars
            mycars$cyl <- factor(mycars$cyl, levels = c("4","6","8"))
            mycars$engine_shape <- factor(c("0" = "V-shaped", "1" = "straight")[as.character(mycars$vs)])
            mycars$transmission <- factor(c("0" = "automatic", "1" = "manual")[as.character(mycars$am)])
            i <- match(name, .info$name)
            if (verbose) {
                message("`", name, "`: ", .info[i, "description"])
            }
            assign(x = "mycars", value = mycars, envir = envir)
            .info[i, ]
        },
        "mycars_summary" = {
            mycars_summary <- datasets::mtcars |>
                dplyr::group_by(.data$cyl) |>
                dplyr::summarise(hp_avg = mean(.data$hp),
                          hp_sem = sd(.data$hp) / sqrt(length(.data$hp)),
                          .groups = "drop")
            i <- match(name, .info$name)
            if (verbose) {
                message("`", name, "`: ", .info[i, "description"])
            }
            assign(x = "mycars_summary", value = mycars_summary, envir = envir)
            .info[i, ]
        },
        "2023-05" = {
            df <- rbind(exampleData("mycars", envir, verbose),
                        exampleData("mycars_summary", envir, verbose))
            df
        }
    )

    invisible(df)

}
