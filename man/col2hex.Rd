% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/col2hex.R
\name{col2hex}
\alias{col2hex}
\title{Get hex color}
\usage{
col2hex(col, alpha = 255)
}
\arguments{
\item{col}{a \code{character}, \code{integer} or vector of both types containing the names of the colors or colors as integers.}

\item{alpha}{a numerical value in the range [0,1] or [0,255] that indicates the 
transparency of the color(s). If the given values are between 0 and 1, they 
are mapped to be between 0 and 255. An alpha value of 1 assumes the [0,1] range
and provides maximum color. The default is set to 255.}
}
\value{
a \code{character} or \code{character} vector with the hex colors.
}
\description{
The function returns a color in hex form given a valid name of a color in R.
}
\examples{
y <- rnorm(1000,0,1)
cols <- rep("red", length(y))
alpha <- seq(0,1,length.out=length(y))
hexcols <- col2hex(cols, alpha)
plot(1:length(y), y, bg=hexcols, pch=21)

y <- rnorm(1000,0,1)
cols <- rep("red", length(y))
alpha <- seq(0,255,length.out=length(y))
hexcols <- col2hex(cols, alpha)
plot(1:length(y), y, bg=hexcols, pch=21)

}
\author{
Dania Machlab
}
