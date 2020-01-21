#' @title Check for valid color specification.
#'
#' @description Check if the elements of a \code{numeric} or \code{character} vector
#'     are valid colors understood by R.
#'
#' @author Michael Stadler
#'
#' @param x \code{numeric} or \code{character} vector with color specifications.
#'
#' @details This is an internal function, used for example by \code{valueToColor}.
#'     The elements of \code{x} are check for valid color specifications by
#'     passing them to \code{col2rgb} and checking for an exception.
#'
#' @return A \code{logical} vector of the same length as \code{x}.
#'
#' @examples
#' .isValidColor(c("black", "blackk"))
#' .isValidColor(c(-1, 1))
#'
#' @keywords internal
#' 
#' @importFrom grDevices col2rgb
.isValidColor <- function(x)
    vapply(X = x, FUN = function(y) tryCatch(is.matrix(grDevices::col2rgb(y)),
                                             error = function(e) FALSE),
           FUN.VALUE = logical(1), USE.NAMES = FALSE)


#' @title Map numerical values to colors.
#'
#' @description \code{valueToColor} takes a numerical vector and maps each value
#'     to an R color string.
#'
#' @author Michael Stadler
#'
#' @param x \code{numeric} vector with values to be mapped to colors.
#' @param rng \code{numeric(2)} giving the range of values to be mapped to colors.
#'     By default, this will be the range of finite values in \code{x}.
#' @param col vector with R colors defining the palette (must be a valid argument
#'     to \code{\link[grDevices]{col2rgb}}.
#' @param NA.col Single R color to use for \code{NA} values in \code{x}.
#' @param alpha \code{NULL} (default) or \code{numeric(1)} between 0 and 255,
#'     giving the alpha channel value for the colors (0 = fully transparent, 255 = fully opaque).
#'     \code{NULL} will use fully opaque colors (\code{alpha = 255}). \code{alpha} is
#'     ignored if \code{col} already contain colors with defined alpha values.
#'
#' @details The values in \code{[rng[1], rng[2]]} will be linearly mapped to the
#'     color palette defined by \code{col}. Any values in \code{x} less (greater)
#'     than \code{rng[1]} (\code{rng[2]}) will be assigned the same color as
#'     \code{rng[1]} (\code{rng[2]}).
#'
#' @return A \code{character} vector of the same length of \code{x} with
#'     R colors in hexadecimal string-encoded RGB format.
#'
#' @seealso \code{\link[grDevices]{colorRamp}} and \code{\link[grDevices]{rgb}}
#'     for the functions called by \code{valueToColor}.
#'
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' cols <- valueToColor(x + y)
#' plot(x, y, pch = 20, col = cols, main = "default")
#'
#' @importFrom grDevices col2rgb colorRamp rgb
#'
#' @export
valueToColor <- function(x, rng = range(x, na.rm = TRUE),
                         col = c("#5E4FA2","#3288BD","#66C2A5","#ABDDA4","#E6F598",
                                 "#FFFFBF","#FEE08B","#FDAE61", "#F46D43","#D53E4F",
                                 "#9E0142"),
                         NA.col = "lightgray", alpha = NULL) {
    # arguments
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(rng) && length(rng) == 2L && rng[1] < rng[2])
    stopifnot(all(.isValidColor(col)))
    stopifnot(.isValidColor(NA.col) && length(NA.col) == 1L)
    stopifnot(is.null(alpha) ||
              (is.numeric(alpha) && length(alpha) == 1L && alpha >= 0 && alpha <= 255))

    # add alpha channel where needed
    colMatrix <- grDevices::col2rgb(col = col, alpha = TRUE)
    if (!is.null(alpha)) {
        if (all(colMatrix[4,] == 255))
            colMatrix[4,] <- as.integer(alpha)
        else
            warning("ignoring 'alpha', as 'col' already specifies alpha channels")
    }
    col <- grDevices::rgb(colMatrix[1,], colMatrix[2,], colMatrix[3,],
                          colMatrix[4,], maxColorValue = 255)

    # truncate x to rng
    i <- !is.na(x)
    x[i & x < rng[1]] <- rng[1]
    x[i & x > rng[2]] <- rng[2]

    # map x to colors
    colfunc    <- grDevices::colorRamp(col, alpha = TRUE)
    xcolMatrix <- colfunc((x[i] - rng[1]) / (rng[2] - rng[1]))
    xcol       <- rep(NA.col, length(x))
    xcol[i]    <- grDevices::rgb(xcolMatrix[,1], xcolMatrix[,2], xcolMatrix[,3],
                                 xcolMatrix[,4], maxColorValue = 255)

    return(xcol)
}
