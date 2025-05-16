#' Log ratio estimator with bootstrap confidence intervals
#'
#' Estimates the geometric mean of the log-ratio of two paired samples with 
#' a bias correction and computes a bootstrap confidence interval.
#'
#' @param x A numeric vector of positive values (e.g., treatment group).
#' @param y A numeric vector of positive values (e.g., control group). Must be the same length as `x`.
#' @param conf.level Confidence level for the bootstrap confidence interval. Defaults to 0.95.
#' @param n_boot Number of bootstrap samples to draw. Defaults to 1000.
#' @param return_boot Logical; if `TRUE`, the function also returns the full bootstrap distribution (on the ratio scale). Defaults to `FALSE`.
#' @return A list with the following components:
#' \describe{
#'   \item{mu}{Bias-corrected estimate of the geometric mean ratio.}
#'   \item{lci_mu}{Lower bound of the bootstrap confidence interval.}
#'   \item{uci_mu}{Upper bound of the bootstrap confidence interval.}
#'   \item{log_mu}{Mean of the log-ratios (i.e., log(x) - log(y)).}
#'   \item{n}{Sample size.}
#'   \item{bootstrap_distribution}{(Optional) Matrix of bootstrap estimates on the ratio scale, returned if `return_boot = TRUE`.}
#' }
#' @export

log_ratio_estimator_bootstrap <- function(x, y, conf.level = 0.95, n_boot = 1000, return_boot = FALSE) {
  
  n <- length(x)
  log_mu <- log(x) - log(y)
  mean_log <- mean(log_mu)
  var_log <- var(log_mu)
  
  # Bias-corrected geometric mean ratio
  mu <- exp(mean_log + 0.5 * var_log / n)
  
  mean_var <- function(x) {
    cbind(mean(x), var(x))
  }
  
  # Bootstrap sampling
  boot_means <- replicate(n_boot, {
    indices <- sample(1:n, size = n, replace = TRUE)
    mean_var(log(x[indices]) - log(y[indices]))
  }, 
  simplify = TRUE)
  
  boot_means <- t(boot_means)
  
  # Calculate bias correction factor for each sample
  boot_means <- cbind(boot_means, 0.5 * boot_means[, 2] / n)
  
  # Bootstrap confidence interval on the ratio scale
  ci_bounds <- quantile(
    exp(boot_means[, 1] + boot_means[, 3]), 
    probs = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  )
  
  result <- list(
    mu = mu,
    lci_mu = ci_bounds[1],
    uci_mu = ci_bounds[2],
    log_mu = mean_log,
    n = n
  )
  
  if(return_boot) {
    result$bootstrap_distribution <- exp(boot_means)
    result$bootstrap_distribution <- cbind(result$bootstrap_distribution, boot_means)
    
    names(result$bootstrap_distribution) <- c("mean", "var", "bc_factor", "bc_mean")
  }
  
  return(result)
}
