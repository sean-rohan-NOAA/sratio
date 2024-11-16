#' Predict Count Based on Selectivity Ratio
#'
#' Calculate the predicted count for one group based on the observed count of the other group, adjusted by the selectivity ratio (`s12`), sampling factor, and effort. If `est_count` is `"count1"`, it estimates the value of `count1` based on `count2` and vice versa.
#'
#' @param est_count Character indicating whether to estimate `"count1"` or `"count2"`.
#' @param count1 Numeric vector of observed counts of size/age in group 1.
#' @param count2 Numeric vector of observed counts of size/age in group 2.
#' @param sampling_factor1 Numeric vector of sampling factors for group 1 (default = 1).
#' @param sampling_factor2 Sampling factor for group 2 (default = 1).
#' @param effort1 Effort for group 1 (default = 1).
#' @param effort2 Effort for group 2 (default = 1).
#' @param s12 Numeric vector of selectivity ratios between group 1 and group 2, where s12 = selectivity1/selectivity2. Values of `0` are automatically adjusted to a small number (`1e-12`) to avoid dividing by zero.
#' @return  Predicted count for the specified group (`count1` or `count2`), rounded to six decimal places.
#' @export

sratio_predict_count <- function(est_count, count1, count2, sampling_factor1 = 1, sampling_factor2 = 1, effort1 = 1, effort2 = 1, s12) {
  
  s12[s12 == 0] <- 1e-12
  
  if(est_count == "count2") {
    count_fit <- (count1 * sampling_factor1 / effort1) / s12 * effort2 / sampling_factor2
  }
  
  if(est_count == "count1") {
    count_fit <- (count2 * sampling_factor2 / effort2) * s12 * effort1 / sampling_factor1
  }
  
  count_fit <- round(count_fit, digits = 6)
  
  return(count_fit)
  
}