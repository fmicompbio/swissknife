#' @title Map numerical values to colors.
#'
#' @description \code{colByValue} takes a numerical vector and maps each value
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
#' @param alpha \code{numeric(1)} between 0 and 255, giving the alpha channel
#'     value for the colors (0 = fully transparent, 255 = fully opaque,
#'     ignored if \code{col} already contain alpha values).
#'
#' @details The values in \code{[rng[1], rng[2]]} will be linearly mapped to the
#'     color palette defined by \code{col}. Any values in \code{x} less (greater)
#'     than \code{rng[1]} (\code{rng[2]}) will be assigned the same color as
#'     \code{rng[1]} (\code{rng[2]}).
#'
#' @return A \code{character} vector of the same length of \code{x} with
#'     R colors in hexadecimal string-encoded RGB format.
#'
#' @seealso \code{\link[grDevices]{colRamp}} and \code{\link[grDevices]{rgb}}
#'     for the functions called by \code{colByValue}.
#'
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' cols <- colByValue(x + y)
#' plot(x, y, pch = 20, col = cols, main = "default")
#'
#' @importFrom grDevices colRamp rgb
#'
#' @export
colByValue <- function(x, rng = range(x, na.rm = TRUE),
                       col = c("#5E4FA2","#3288BD","#66C2A5","#ABDDA4","#E6F598",
                               "#FFFFBF","#FEE08B","#FDAE61", "#F46D43","#D53E4F",
                               "#9E0142"),
                       NA.col = "lightgray", alpha = 255) {
    # helper function
    isValidColor <- function(y)
        vapply(X = y, FUN = function(Y) tryCatch(is.matrix(col2rgb(Y)),
                                                 error = function(e) FALSE),
               FUN.VALUE = logical(1), USE.NAMES = FALSE)

    # arguments
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(rng) && length(rng) == 2L && rng[1] < rng[2])
    stopifnot(all(isValidColor(col)))
    stopifnot(isValidColor(NA.col) && length(NA.col) == 1L)
    stopifnot(is.numeric(alpha) && length(alpha) == 1L && alpha >= 0 && alpha <= 255)

    # add alpha channel where needed
    r <- grDevices::col2rgb(col = col, alpha = TRUE)
    if (all(r[4,] == 255))
        r[4,] <- as.integer(alpha)
    else if (alpha != 255)
        warning("ignoring 'alpha', as 'col' already contains alpha channels")
    col <- grDevices::rgb(r[1,], r[2,], r[3,], r[4,], maxColorValue = 255)

    # truncate x to rng
    i <- !is.na(x)
    x[i & x < rng[1]] <- rng[1]
    x[i & x > rng[2]] <- rng[2]

    # map x to colors
    colfunc <- grDevices::colorRamp(col, alpha = TRUE)
    r <- colfunc((x[i] - rng[1]) / (rng[2] - rng[1]))
    rr <- rep(NA.col, length(x))
    rr[i] <- grDevices::rgb(r[,1], r[,2], r[,3], r[,4], maxColorValue = 255)

    return(rr)
}
