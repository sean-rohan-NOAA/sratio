#' Inverse logit
#' 
#' Inverse logit
#' 
#' @param x Numeric vector
#' @export

inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}
