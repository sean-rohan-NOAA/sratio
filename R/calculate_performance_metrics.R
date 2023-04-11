#' Calculate RMSLE, MRE, MAE, R Sq., Bias
#' 
#' Compare series x to y using performance metrics from Seegers et al. (2018). Root mean square log error, mean relative error, mean absolute error, R-squared, and mean bias.
#' 
#' @param x numeric vector for the first series
#' @param y numeric vector for the second series
#' @param add_for_log Value to add to x and y that require log transformations.
#' @export


calculate_performance_metrics <- function(x, y, add_for_log = 1) {
  
  root_mean_square_log_error <- mean(sqrt((log(x + add_for_log) - log(y + add_for_log))^2))
  mean_relative_error <- mean(abs(x[y != 0 & x != 0] - y[y != 0 & x != 0]) / y[y != 0 & x != 0])
  mean_absolute_error <- 10^mean(abs(log10(x + add_for_log) - log10(y + add_for_log)))
  rsq <- cor(y, x)^2
  mean_bias <- 10^mean(((log10(x + add_for_log) - log10(y + add_for_log))))
  
  return(list(bias = mean_bias,
              rmsle = root_mean_square_log_error,
              mre = mean_relative_error,
              mae = mean_absolute_error,
              rsq = rsq))
  
}