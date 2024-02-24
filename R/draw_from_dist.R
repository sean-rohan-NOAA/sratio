#' Draw samples from a distribution
#' 
#' @param x A list containing arguments passed to the distribution of choice.
#' @details
#' x should be a list containing an argument named distribution denoting the distribution to use (e.g. 'normal') and parameters passed to the distributional function for the distribution (e.g. the normal distribution specifying the use of rnorm should have mean and sd arguments)
#' @export

draw_from_dist <- function(x) {
  
  stopifnot("draw_from_dist: x must be a list" = class(x) == "list")
  stopifnot("draw_from_dist: x must include a distribution argument denoting the distrubution to use." = "distribution" %in% names(x))
  
  fn <- c("rnorm")[match(x$distribution, table = "normal")]
  
  if(!("n" %in% names(x))) {
    x$n <- 1
  }
  
  output <- do.call(what = fn, args = x[-which(names(x) == "distribution")])
  
  return(output)
  
}