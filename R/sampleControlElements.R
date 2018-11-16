#' @title Sample control elements that match a target distribution.
#'
#' @description Randomly sample from a set of control (background) elements,
#'     such that the selected elements are similarly distributed as a given
#'     set of target (foreground) elements.
#'
#' @author Michael Stadler
#'
#' @param x `numeric` vector (or list of `numeric` vectors). `idxTarget` and `idxControl`
#'     refer to the elements of `x`. If `x` is a list, all elements must have the same length.
#' @param idxTarget `numeric` or `logical` vector specifying the elements in `x`
#'     that define the target distribution to be matched by the control elements.
#' @param idxControl `numeric` or `logical` vector specifying the complete set of
#'     possible control elements in `x` (default: all that are not in `idxTarget`),
#'     from which a subset is to be sampled.
#' @param nbins `numeric(1)` or `numeric(length(x))` if `x` is a list, specifying the
#'     number of bins to group the values of x into. Higher numbers of bins will
#'     increase the match to the target distribution(s), but may fail if there are
#'     few elements to sample from.
#' @param oversample The number of control elements to sample for each target element.
#'
#' @return `numeric` vector with `round(length(idxTarget) * oversample)` elements,
#'     specifying the index (positions) of the sampled control elements.
#'
#' @examples
#' x <- c(runif(1000, min = 0, max = 10),
#'        rnorm(200, mean = 5, sd = 1))
#' s <- sampleControlElements(x, idxTarget = 1001:1200, idxControl = 1:1000)
#' par(mfrow=c(2,2))
#' h <- hist(x, breaks = 20, main = "all")
#' hist(x[1:1000], breaks = h$breaks, main = "all control")
#' hist(x[1001:1200], breaks = h$breaks, main = "target")
#' hist(x[s], breaks = h$breaks, main = "sampled control")
#'
#' @export
sampleControlElements <- function(x, idxTarget, idxControl = NULL, nbins = 50, oversample = 1.0) {
    # ... check arguments
    if (!is.list(x)) {
        if (is.numeric(x)) {
            x <- list(x)
        } else {
            stop("`x` needs to be either a `list` or a `numeric` vector")
        }
    } else if (!all(unlist(lapply(x, is.numeric)))) {
        stop("all elements of the list `x` have to be `numeric`")
    }
    stopifnot(is.list(x) && all(lengths(x) == length(x[[1]])))
    if (is.logical(idxTarget))
        idxTarget <- which(idxTarget)
    if (is.null(idxControl))
        idxControl <- seq_along(x[[1]])[-idxTarget]
    if (is.logical(idxControl))
        idxControl <- which(idxControl)
    stopifnot(is.numeric(idxTarget) && is.numeric(idxControl))
    N <- length(x)
    n <- length(x[[1]])
    if (length(nbins) == 1L && N > 1)
        nbins <- rep(nbins, N)
    stopifnot(is.numeric(nbins) && length(nbins) == N)
    stopifnot(min(idxTarget,idxControl) >= 1 && max(idxTarget,idxControl) <= n)

    # ... bin target (convert to factor to retain unobserved levels)
    binsL <- lapply(seq_along(x), function(i) seq(min(x[[i]]), max(x[[i]]), length.out = nbins[i] + 1L))
    target.binsL <- lapply(seq_along(x), function(i) factor(findInterval(x[[i]][idxTarget], binsL[[i]], all.inside = TRUE),
                                                            levels = seq.int(nbins[i])))
    target.bins <- table(target.binsL)

    # ... bin control
    control.binsL <- lapply(seq_along(x), function(i) factor(findInterval(x[[i]][idxControl], binsL[[i]], all.inside = TRUE),
                                                             levels = seq.int(nbins[i])))
    control.bins <- table(control.binsL)
    control.bins.tab <- do.call(cbind, control.binsL)

    # ... select control elements using target-based probabilities
    control.prob <- target.bins[control.bins.tab] / control.bins[control.bins.tab]
    sel <- sample(x = idxControl, size = round(length(idxTarget)*oversample),
                  prob = control.prob)
    return(sel)
}
