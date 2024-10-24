#' Non-Parametric Bootstrap for Estimating Mean Bias
#'
#' This function performs a non-parametric bootstrap to estimate the bias between two datasets, \code{x1} and \code{x2}, by computing a mean bias across multiple resamples. The user can specify the number of bootstrap samples, whether to resample with replacement, the scale of the calculation, and an optional constant to add to the data values before bias estimation.
#'
#' @param x1 A numeric vector. First dataset.
#' @param x2 A numeric vector. Second dataset, must be the same length as \code{x1}.
#' @param n_samples An integer. Number of bootstrap samples to generate.
#' @param add_constant A numeric value. A constant added to both \code{x1} and \code{x2} before computing the bias. Default is \code{1}.
#' @param scale A character string. The scale for bias calculation, either \code{"log10"} for logarithmic (base 10) scale or \code{"none"} for linear scale. Default is \code{"log10"}.
#' @param seed An optional integer. Random seed for reproducibility. Default is \code{NULL}, meaning no specific seed is set.
#'
#' @details The bias is computed as:
#' \itemize{
#'   \item When \code{scale = "log10"}: \eqn{10^{\text{mean}(\log_{10}(x1 + \text{add_constant}) - \log_{10}(x2 + \text{add_constant}))}}
#'   \item When \code{scale = "none"}: \eqn{\text{mean}((x1 + \text{add_constant}) - (x2 + \text{add_constant}))}
#' }
#'
#' The function resamples \code{x1} and \code{x2} \code{n_samples} times, computes the bias for each sample, and returns a vector of the bias values.
#'
#' @return A numeric vector of the computed mean bias for each bootstrap sample.
#'
#' @examples
#' # Example usage:
#' x1 <- c(2, 3, 5, 7)
#' x2 <- c(1, 2, 3, 4)
#' bootstrap_mean_bias(x1, x2, n_samples = 1000, add_constant = 1, scale = "log10", seed = 42)
#'
#' @export

bootstrap_mean_bias <- function(x1, x2, n_samples, add_constant = 1, scale = "log10", seed = NULL) {
  
  stopifnot("bootstrap_mean_bias: x1 and x2 must be the same length" = length(x1) == length(x2))
  
  mean_bias <- numeric()
  
  if(scale == "log10") {
    bias_fn <- function(a, b, c) {
      return(10^mean(log10(a + c) - log10(b + c)))
    }
  }
  
  if(scale == "none") {
    bias_fn <- function(a, b, c) {
      return(mean((a + c) - (b + c)))
    }
  }
  
  set.seed(seed)
  
  for(ii in 1:n_samples) {
    ind <- sample(x = 1:length(x1), size = length(x1), replace = TRUE)
    sample_bias <- bias_fn(a = x1[ind], b = x2[ind], c = add_constant)
    mean_bias <- c(mean_bias, sample_bias)
  }
  
  return(mean_bias)
  
}