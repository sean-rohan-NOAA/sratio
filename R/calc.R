#' Calculate Bias
#'
#' This function calculates the bias between estimated and observed values.
#'
#' @param est A numeric vector of estimated values.
#' @param obs A numeric vector of observed values.
#' @param const A numeric constant added to the estimated and observed values before calculation. Default is 0.
#' @param type A function to summarize the bias, such as \code{mean} or \code{median}. Default is \code{mean}.
#' @return A numeric value representing the bias.
#' @examples
#' est <- c(1, 2, 3)
#' obs <- c(1.1, 2.1, 2.9)
#' calc_bias(est, obs)
#' @export
calc_bias <- function(est, obs, const = 0, type = mean) {
  10^type(log10(est + const) - log10(obs + const))
}

#' Calculate Root Mean Square Error (RMSE)
#'
#' This function calculates the root mean square error between estimated and observed values.
#'
#' @param est A numeric vector of estimated values.
#' @param obs A numeric vector of observed values.
#' @param type A function to summarize the error, such as \code{mean} or \code{median}. Default is \code{mean}.
#' @return A numeric value representing the RMSE.
#' @examples
#' est <- c(1, 2, 3)
#' obs <- c(1.1, 2.1, 2.9)
#' calc_rmse(est, obs)
#' @export
calc_rmse <- function(est, obs, type = mean) {
  sqrt(type((est - obs)^2))
}

#' Calculate Mean Absolute Error (MAE)
#'
#' This function calculates the mean absolute error between estimated and observed values.
#'
#' @param est A numeric vector of estimated values.
#' @param obs A numeric vector of observed values.
#' @param const A numeric constant added to the estimated and observed values before calculation. Default is 0.
#' @param type A function to summarize the error, such as \code{mean} or \code{median}. Default is \code{mean}.
#' @return A numeric value representing the MAE.
#' @examples
#' est <- c(1, 2, 3)
#' obs <- c(1.1, 2.1, 2.9)
#' calc_mae(est, obs)
#' @export
calc_mae <- function(est, obs, const = 0, type = mean) {
  10^type(abs(log10(est + const) - log10(obs + const)))
}

#' Calculate Mean Relative Error (MRE)
#'
#' This function calculates the mean relative error between estimated and observed values.
#'
#' @param est A numeric vector of estimated values.
#' @param obs A numeric vector of observed values.
#' @param type A function to summarize the error, such as \code{mean} or \code{median}. Default is \code{mean}.
#' @return A numeric value representing the MRE.
#' @examples
#' est <- c(1, 2, 3)
#' obs <- c(1.1, 2.1, 2.9)
#' calc_mre(est, obs)
#' @export
calc_mre <- function(est, obs, type = mean) {
  re <- (est - obs) / obs * 100
  length_inf <- sum(is.infinite(re))
  if (length_inf > 0) {
    warning("calc_mre: Ignoring ", length_inf, " observations where obs = 0.")
  }
  type(re[!is.infinite(re)], na.rm = TRUE)
}
