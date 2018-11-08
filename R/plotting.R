# functions related to plotting / visualization

#' @title Create a bitmap-rendered plot.
#'
#' @description `plotBitScatter` is a wrapper around `plot` which renders the
#'     plot area as a bitmap (png), but keeps all other elements (axes, labels, etc.)
#'     as vector elements. This is especially useful for keeping the size of PDF files
#'     with scatter plots with many elements small, while retaining editability of axes.
#'
#' @author Michael Stadler
#'
#' @param x `numeric` vector with x-coordinates of points.
#' @param y `numeric` vector with y-coordinates of points (same length as `x`).
#' @param ... any further arguments to be passed to `plot`
#' @param densCols `logical(1)`. If `TRUE` and `col` is not given as an additional
#'     argument, then the local density of points will be used as colors, using the
#'     palette spanned by the colors in `colpal`.
#' @param colpal vector of colors defining the palette for automatic density-based coloring.
#' @param xpixels the number of pixels in the x dimension used for rendering
#'     the plotting area. The number of pixels in the y dimension are calculated
#'     as `xpixels * par('pin')[2] / par('pin')[1]`, such that the aspect ratio of
#'     the current plotting region is observed.
#'
#' @details `xpixels` controls the resolution of the rendered plotting area. In order
#'     to keep circular plotting symbols circlular (e.g. `pch = 1`), `ypixels` is
#'     automatically calculated by adjusting `xpixels` to the aspect ratio of the
#'     current plotting area. If the plotting device is rescaled after calling
#'     `plotBitScatter`, circular plotting symbols may become skewed.
#'
#' @return `NULL` (invisibly)
#'
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' par(mfrow=c(1,2))
#' plotBitScatter(x, y, main = "bitmap")
#' plot(x, y, main = "default")
#'
#' @importFrom grDevices png as.raster densCols colorRampPalette dev.off
#' @importFrom png readPNG
#' @importFrom graphics rasterImage par grconvertX plot box
#' @importFrom KernSmooth bkde2D
#'
#' @export
plotBitScatter <- function(x, y, ..., densCols=TRUE,
                           colpal=c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),
                           xpixels=1000) {
    # digest arguments
    args <- list(x = x, y = y, ...)
    if (!"xlab" %in% names(args))
        args[["xlab"]] <- "x"
    if (!"ylab" %in% names(args))
        args[["ylab"]] <- "y"
    if (!"col" %in% names(args) && densCols)
        args[["col"]] <- grDevices::densCols(x, y, nbin = 512, colramp = grDevices::colorRampPalette(colpal))

    # get aspect ratio of current plotting region (will open one if none is open yet)
    ar <- graphics::par("pin")[2] / graphics::par("pin")[1]

    # adjust default png device pointsize (12) to xpixels
    ps <- 12 / graphics::grconvertX(par("pin")[1], from = "inches", to = "device") * xpixels

    # create scatter plot without annotation into a temporary png file
    # REMARK: use grid::grid.cap() and on-screen device instead of temporary png file?
    tf <- tempfile(fileext = ".png")
    on.exit(unlink(tf))
    grDevices::png(tf, width = xpixels, height = xpixels * ar, pointsize = ps)
    graphics::par(mar = c(0,0,0,0))
    do.call("plot", c(list(ann = FALSE, axes = FALSE), args[!(names(args) %in% c("ann","axes"))]))
    grDevices::dev.off()

    # read in png and create raster object
    im <- grDevices::as.raster(png::readPNG(tf))

    # create main plot and add raster image to it
    do.call(graphics::plot, c(list(type = "n"), args[!(names(args) %in% c("type"))]))
    usr <- graphics::par("usr")
    graphics::rasterImage(image = im, xleft = usr[1], ybottom = usr[3],
                          xright = usr[2], ytop = usr[4], interpolate = TRUE)
    graphics::box()
    return(invisible(NULL))
}
