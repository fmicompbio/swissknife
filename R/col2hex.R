#' @title Get hex color
#'
#' @author Dania Machlab
#' 
#' @description The function returns a color in hex form given a valid name of a color in R.
#' 
#' @param col a \code{character}, \code{integer} or vector of both types containing the names of the colors or colors as integers.
#' @param alpha a numerical value in the range [0,1] or [0,255] that indicates the 
#'   transparency of the color(s). If the given values are between 0 and 1, they 
#'   are mapped to be between 0 and 255. An alpha value of 1 assumes the [0,1] range
#'   and provides maximum color. The default is set to 255.
#'   
#' @return a \code{character} or \code{character} vector with the hex colors. 
#' 
#' @examples 
#' y <- rnorm(1000,0,1)
#' cols <- rep("red", length(y))
#' alpha <- seq(0,1,length.out=length(y))
#' hexcols <- col2hex(cols, alpha)
#' plot(1:length(y), y, bg=hexcols, pch=21)
#' 
#' y <- rnorm(1000,0,1)
#' cols <- rep("red", length(y))
#' alpha <- seq(0,255,length.out=length(y))
#' hexcols <- col2hex(cols, alpha)
#' plot(1:length(y), y, bg=hexcols, pch=21)
#' 
#' @importFrom grDevices col2rgb rgb
#'
#' @export
col2hex <- function(col, alpha=255) {
  
  # checks
  if(any(!swissknife:::.isValidColor(col))) {
    stop("the color provided is not a valid color in R")
  }
  if(!is.numeric(alpha)){
    stop("alpha must be numeric")
  }
  if(length(alpha)!=1 & length(alpha)!=length(col)){
    stop("alpha must be the same length as the col vector if the alphas ought to be specific to the elements in col")
  }
  
  # If alpha is in [0,1], map it to [0,255]
  if(all(alpha>=0 & alpha<=1)){
    alpha <- alpha*255
  }
  
  # check alpha range
  if(any(!(alpha>=0 & alpha<=255))){
    stop("alpha must be in the range [0,1] or [0,255]")
  }
  
  # get rgb
  rgb <- grDevices::col2rgb(col)
  
  # get and return hex
  grDevices::rgb(red = rgb["red", ], green = rgb["green", ], blue = rgb["blue", ], maxColorValue=255, alpha = alpha)
  
}




