#' Log ratio estimator with bootstrap confidence intervals
#'
#' Estimates the geometric mean of the log-ratio of two paired samples
#'
#' @param x A numeric vector of positive values.
#' @param y A numeric vector of positive values.
#' @param conf_level Confidence level for the bootstrap confidence interval. Defaults to 0.95.
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

log_ratio_estimator_bootstrap <- function(x, y, conf_level = 0.95, n_boot = 1000, return_boot = FALSE) {
  
  # Original log-ratio estimate
  log_ratios <- log(x / y)
  theta_hat <- mean(log_ratios)
  
  # Bootstrap sampling
  n <- length(x)
  boot_estimates <- numeric(n_boot)
  
  for (b in 1:n_boot) {
    indices <- sample(1:n, size = n, replace = TRUE)
    boot_log_ratios <- log(x[indices] / y[indices])
    boot_estimates[b] <- mean(boot_log_ratios)
  }
  
  # Standard error and CI
  se <- sd(boot_estimates)
  alpha <- 1 - conf_level
  ci <- quantile(boot_estimates, probs = c(alpha/2, 1 - alpha/2))
  
  # Return results
  return(list(
    estimate = theta_hat,
    bootstrap_se = se,
    conf_interval = ci,
    boot_estimates = boot_estimates
  ))
  
  return(result)
}
