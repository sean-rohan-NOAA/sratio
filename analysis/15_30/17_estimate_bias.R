# Calculate CPUE ratio, bias, mean absolute error, and root mean square error 

library(sratio)

seed <- 1729730909

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

# Load built-in data sets
catch_df <- sratio::data_1530$catch #|>
  # dplyr::filter(CRUISE %in% use_cruises)

haul_df <- sratio::data_1530$haul #|>
  # dplyr::filter(CRUISE %in% use_cruises)

sp_codes <- sort(unique(catch_df$SPECIES_CODE))

# bias_table <- catch_df |>
#   dplyr::inner_join(haul_df) |>
#   dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
#   dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
#   tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
#   dplyr::inner_join(readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
#                       dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
#                       unique() |>
#                       dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
#   dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
#                 CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
#   dplyr::mutate(ratio = CPUE_15/CPUE_30,
#                 log_error = log10(CPUE_30+1)-log10(CPUE_15+1),
#                 abs_error = abs(log10(CPUE_30+1)-log10(CPUE_15+1)),
#                 sq_error = (CPUE_30-CPUE_15)^2) |>
#   dplyr::group_by(SPECIES_CODE) |>
#   dplyr::reframe(MEAN_RATIO = mean(ratio),
#                    BIAS = 10^(mean(log_error)),
#                    MAE = 10^mean(abs_error),
#                    RMSE = sqrt(mean(sq_error))) |>
#   dplyr::mutate(COMMON_NAME = sratio:::species_code_label(SPECIES_CODE, type = "common_name")) |>
#   dplyr::select(COMMON_NAME, BIAS, MEAN_RATIO, MAE, RMSE)

bias_samples <- catch_df |>
  dplyr::inner_join(haul_df) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
                CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::reframe(BIAS = bootstrap_mean_bias(CPUE_30, CPUE_15, n_samples = 1000, replace = TRUE, add_constant = 1, scale = "log10", seed = seed)) |>
  dplyr::mutate(COMMON_NAME = sratio:::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::select(COMMON_NAME, BIAS)

bias_quantiles <- bias_samples |>
  dplyr::group_by(COMMON_NAME) |>
  dplyr::summarise(q025 = quantile(BIAS, 0.025),
                   q250 = quantile(BIAS, 0.25),
                   q500 = quantile(BIAS, 0.5),
                   q750 = quantile(BIAS, 0.75),
                   q975 = quantile(BIAS, 0.975),
                   type = "Bias")


lines <- c("CPUE comparison between 15 and 30 minute tows\n", "BIAS > 1 = 30 minutes higher\n\n\n")
cat(lines, file = here::here("analysis", "15_30", "plots", "total_cpue_fit", "bias_table.csv"))

write.table(bias_table, 
            file = here::here("analysis", "15_30", "plots", "total_cpue_fit", "bias_table.csv"), 
            append = TRUE, 
            row.names = FALSE, 
            sep = ",")
