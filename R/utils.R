#' Connect to Oracle
#'
#' Connect to set up a connection to Oracle using RODBC and getPass
#'
#' @param schema Connection. Default = 'AFSC'
#' @export

get_connected <- function(schema='AFSC'){(echo=FALSE)
  username <- getPass::getPass(msg = "Enter your ORACLE Username: ")
  password <- getPass::getPass(msg = "Enter your ORACLE Password: ")
  channel  <- RODBC::odbcConnect(paste(schema),paste(username),paste(password), believeNRows=FALSE)
}


#' Calculate weighted mean
#'
#' Calculate the weighted mean of a numeric vector, given corresponding weights.
#' 
#' @param x A numeric vector of values for which the weighted mean is to be calculated.
#' @param w A numeric vector of weights corresponding to each value.
#' @param na_rm Should NA's be ignored?
#'
#' @return The weighted mean of the values.
#' 
#' @examples
#' x <- c(1, 2, 3, 4)
#' w <- c(0.1, 0.2, 0.3, 0.4)
#' weighted_mean(x, w)
#' @export

weighted_mean <- function(x, w, na_rm = FALSE) {
  
  stopifnot("weighted_mean: The length of x and w must be the same." = length(x) == length(w))
  
  weighted_sum <- sum(x * w, na.rm = na_rm)
  total_weight <- sum(w, na.rm = na_rm)
  
  if(total_weight == 0) {
    stop("The sum of w must not be zero.")
  }
  
  return(weighted_sum / total_weight)
}