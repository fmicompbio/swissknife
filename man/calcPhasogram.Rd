% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phasograms.R
\name{calcPhasogram}
\alias{calcPhasogram}
\title{Calculate phasograms (same strand alignment distances).}
\usage{
calcPhasogram(fname, regions = NULL, rmdup = TRUE, dmax = 3000L)
}
\arguments{
\item{fname}{\code{character} vector with one or several bam files. If
multiple files are given, distance counts from all will be summed.}

\item{regions}{\code{GRanges} object. Only alignments falling into these
regions will be used. If \code{NULL} (the default), all alignments are
used.}

\item{rmdup}{\code{logical(1)} indicating if duplicates should be removed. If
\code{TRUE} (the default), only one of several alignments starting at the
same coordinate is used.}

\item{dmax}{\code{numeric(1)} specifying the maximal distance between same
strand alignments to count.}
}
\value{
\code{integer} vector with \code{dmax} elements, with the element at
  position \code{d} giving the observed number of alignment pairs at that
  distance.
}
\description{
Calculate the frequencies of same strand alignment distances,
  for example from MNase-seq data to estimate nucleosome repeat length.
  Distance calculations are implemented in C++ (\code{\link{calcAndCountDist}})
  for efficiency.
}
\examples{
if (requireNamespace("GenomicAlignments", quietly = TRUE) &&
    requireNamespace("Rsamtools", quietly = TRUE)) {
    bamf <- system.file("extdata", "phasograms", "mnase_mm10.bam",
                        package = "swissknife")
    pg <- calcPhasogram(bamf)
    print(estimateNRL(pg, usePeaks = 1:4)[1:2])
    plotPhasogram(pg, usePeaks = 1:4, xlim = c(0,1000))
}

}
\references{
Phasograms were originally described in Valouev et al., Nature
  2011 (doi:10.1038/nature10002). The implementation here differs in two ways
  from the original algorithms:
  \enumerate{
    \item It does not implement removing of positions that have been seen less
    than \code{n} times (referred to as a \code{n}-pile subset in the paper).
    \item It does allow to retain only alignments that fall into selected
    genomic intervals (\code{regions} argument).
  }
}
\seealso{
\code{\link{estimateNRL}} to estimate the nucleosome repeat length
  from a phasogram, \code{\link{plotPhasogram}} to visualize an annotated
  phasogram, \code{\link{calcAndCountDist}} for low-level distance counting.
}
\author{
Michael Stadler
}
