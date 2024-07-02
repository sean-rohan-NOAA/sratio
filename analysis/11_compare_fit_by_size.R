library(sratio)

sp_codes <- c(21740, 21720, 10110, 10112, 10115, 10130, 10210, 10261, 10285, 471, 68560, 68580, 69322)

central_fun <- median

model_method <- "sccal_model" #"sratio"

# Minimum catch-at-length for 30 minute tows to be used for performance metrics
min_n <- 3

for(ii in 1:length(species_codes)) {
  
  # Load catch-at-length and effort data
  size_30 <- readRDS(here::here("output", "catch_at_length_1530.rds")) |>
    dplyr::filter(TREATMENT == 30, SPECIES_CODE == sp_codes[ii]) |>
    dplyr::ungroup()
  
  size_15 <- readRDS(here::here("output", "catch_at_length_1530.rds")) |>
    dplyr::filter(TREATMENT == 15, SPECIES_CODE == sp_codes[ii]) |>
    dplyr::ungroup()
  
  # Load bootstrap results from best-fitting model
  bootstrap_results_path <- list.files(here::here("output"), 
                                       recursive = TRUE, 
                                       pattern = paste0(model_method, "_bootstrap_results_"), 
                                       full.names = TRUE)
  
  file_index <- grep(pattern = sp_codes[ii], 
                     bootstrap_results_path)
  
  if(length(file_index) < 1) {
    next
  }
  
  # Get mean selectivity ratio
  bootstrap_results <- readRDS(bootstrap_results_path[file_index]) |>
    dplyr::group_by(SIZE_BIN) |>
    dplyr::summarise(s12 = central_fun(s12))
  
  size_15 <- dplyr::select(size_30, MATCHUP, AREA_SWEPT_KM_30 = AREA_SWEPT_KM2) |> 
    unique() |>
    dplyr::inner_join(size_15, by = join_by(MATCHUP)) |>
    dplyr::inner_join(bootstrap_results, by = join_by(SIZE_BIN))
  
  # Predict mean catch-at-length, adjusting for area swept and selectivity
  size_15$PREDICTED_FREQUENCY <- size_15$FREQ_EXPANDED * size_15$AREA_SWEPT_KM_30/size_15$AREA_SWEPT_KM2 * size_15$s12
  
  # Predict mean catch-at-length adjusting for area swept but not selectivity
  size_15$PREDICTED_FREQUENCY_NO_ADJ <- size_15$FREQ_EXPANDED * size_15$AREA_SWEPT_KM_30/size_15$AREA_SWEPT_KM2
  
  compare_fit <- dplyr::select(size_15, MATCHUP, PREDICTED_FREQUENCY, PREDICTED_FREQUENCY_NO_ADJ, 
                               SIZE_BIN, SPECIES_CODE) |>
    dplyr::full_join(size_30, by = join_by(MATCHUP, SIZE_BIN, SPECIES_CODE))  |>
    dplyr::mutate(PREDICTED_FREQUENCY = if_else(is.na(PREDICTED_FREQUENCY), 0, PREDICTED_FREQUENCY),
                  PREDICTED_FREQUENCY_NO_ADJ = if_else(is.na(PREDICTED_FREQUENCY_NO_ADJ), 0, PREDICTED_FREQUENCY_NO_ADJ),
                  FREQ_EXPANDED = if_else(is.na(FREQ_EXPANDED), 0, FREQ_EXPANDED)) |>
    dplyr::filter(FREQ_EXPANDED >= min_n) # Filter records with insufficient samples
  
  compare_fit_long <- tidyr::pivot_longer(compare_fit, 
                                          cols = c("PREDICTED_FREQUENCY", "PREDICTED_FREQUENCY_NO_ADJ"))
  
  p_fit_sratio_adjustment <- ggplot(data = compare_fit,
                                    mapping = aes(x = SIZE_BIN,
                                                  y = PREDICTED_FREQUENCY-FREQ_EXPANDED)) +
    geom_boxplot(mapping = aes(group = SIZE_BIN)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(name = expression("Selectivity ratio "~hat(N[30])-N[30])) +
    scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
    theme_bw()
  
  p_area_only <- ggplot(data = compare_fit,
                        mapping = aes(x = SIZE_BIN,
                                      y = PREDICTED_FREQUENCY_NO_ADJ-FREQ_EXPANDED)) +
    geom_boxplot(mapping = aes(group = SIZE_BIN)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(name = expression("No SR adjustment "~hat(N[30])-N[30])) +
    scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
    theme_bw()
  
  p_fit_vs_obs <- ggplot() +
    geom_point(data = compare_fit,
               mapping = aes(x = PREDICTED_FREQUENCY,
                             y = FREQ_EXPANDED)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_x_continuous(name = expression("Selectivity ratio fit "~hat(N[30]))) + 
    scale_y_continuous(name = expression("Observed "~N[30])) +
    facet_wrap(~SIZE_BIN, scale = "free") +
    theme_bw()
  
  p_bias <- ggplot(data = compare_fit,
                   mapping = aes(x = SIZE_BIN,
                                 y = 10^(log10(PREDICTED_FREQUENCY)-log10(FREQ_EXPANDED)))) +
    geom_boxplot(mapping = aes(group = SIZE_BIN)) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_y_continuous(name = "Bias") +
    scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
    theme_bw()
  
  p_mre <- ggplot(data = compare_fit,
                  mapping = aes(x = SIZE_BIN,
                                y = (PREDICTED_FREQUENCY-FREQ_EXPANDED)/FREQ_EXPANDED*100)) +
    geom_boxplot(mapping = aes(group = SIZE_BIN)) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_y_continuous(name = "MRE (%)") +
    scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
    theme_bw()
  
  # Calculate performance metrics
  dplyr::group_by(compare_fit, SIZE_BIN, SPECIES_CODE) |>
    dplyr::summarise(N_HAUL = n(),
                     MAE = calc_mae(est = PREDICTED_FREQUENCY, 
                                    obs = FREQ_EXPANDED, 
                                    const = 1e-3),
                     MRE = suppressWarnings(calc_mre(est = PREDICTED_FREQUENCY, 
                                    obs = FREQ_EXPANDED)),
                     BIAS = calc_bias(est = PREDICTED_FREQUENCY, 
                                      obs = FREQ_EXPANDED, 
                                      const = 1e-3),
                     R2 = suppressWarnings(cor(PREDICTED_FREQUENCY, 
                              FREQ_EXPANDED, 
                              use = "complete.obs")^2),
                     RMSE = calc_rmse(est = PREDICTED_FREQUENCY, 
                                      obs = FREQ_EXPANDED),
                     MAE_NO_ADJ = calc_mae(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                           obs = FREQ_EXPANDED, 
                                           const = 1e-3),
                     MRE_NO_ADJ = suppressWarnings(calc_mre(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                           obs = FREQ_EXPANDED)),
                     BIAS_NO_ADJ = calc_bias(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = FREQ_EXPANDED, 
                                             const = 1e-3),
                     R2_NO_ADJ = suppressWarnings(cor(PREDICTED_FREQUENCY_NO_ADJ, 
                                     FREQ_EXPANDED, 
                                     use = "complete.obs")^2),
                     RMSE_NO_ADJ = calc_rmse(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = FREQ_EXPANDED),
                     .groups = "keep") |>
    dplyr::mutate(DIFF_RMSE = RMSE - RMSE_NO_ADJ,
                  DIFF_MAE = MAE - MAE_NO_ADJ,
                  DIFF_MRE = abs(MRE) - abs(MRE_NO_ADJ),
                  DIFF_BIAS = abs(BIAS-1) - abs(BIAS_NO_ADJ-1),
                  DIFF_R2 = R2 - R2_NO_ADJ) |>
    as.data.frame()
  
  
  # Estimate total catch in numbers
  
  compare_total_by_haul <- compare_fit |>
    dplyr::group_by(MATCHUP, SPECIES_CODE) |>
    dplyr::summarise(TOTAL_PREDICTED_FREQUENCY = sum(PREDICTED_FREQUENCY),
                     TOTAL_PREDICTED_FREQUENCY_NO_ADJ = sum(PREDICTED_FREQUENCY_NO_ADJ),
                     TOTAL_FREQ_EXPANDED = sum(FREQ_EXPANDED),
                     .groups = "keep")
  
  compare_total <- compare_total_by_haul |> 
    dplyr::group_by(SPECIES_CODE) |>
    dplyr::summarise(N_HAUL = n(),
                     MAE = calc_mae(est = TOTAL_PREDICTED_FREQUENCY, 
                                    obs = TOTAL_FREQ_EXPANDED, 
                                    const = 1e-3),
                     MRE = calc_mre(est = TOTAL_PREDICTED_FREQUENCY, 
                                    obs = TOTAL_FREQ_EXPANDED),
                     BIAS = calc_bias(est = TOTAL_PREDICTED_FREQUENCY, 
                                      obs = TOTAL_FREQ_EXPANDED, 
                                      const = 1e-3),
                     R2 = cor(TOTAL_PREDICTED_FREQUENCY, 
                              TOTAL_FREQ_EXPANDED, 
                              use = "complete.obs")^2,
                     RMSE = calc_rmse(est = TOTAL_PREDICTED_FREQUENCY, 
                                      obs = TOTAL_FREQ_EXPANDED),
                     MAE_NO_ADJ = calc_mae(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                           obs = TOTAL_FREQ_EXPANDED, 
                                           const = 1e-3),
                     MRE_NO_ADJ = calc_mre(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                           obs = TOTAL_FREQ_EXPANDED),
                     BIAS_NO_ADJ = calc_bias(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = TOTAL_FREQ_EXPANDED, 
                                             const = 1e-3),
                     R2_NO_ADJ = cor(TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                     TOTAL_FREQ_EXPANDED, 
                                     use = "complete.obs")^2,
                     RMSE_NO_ADJ = calc_rmse(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = TOTAL_FREQ_EXPANDED),
                     .groups = "keep") |>
    dplyr::mutate(DIFF_RMSE = RMSE - RMSE_NO_ADJ,
                  DIFF_MAE = MAE - MAE_NO_ADJ,
                  DIFF_MRE = abs(MRE) - abs(MRE_NO_ADJ),
                  DIFF_BIAS = abs(BIAS-1) - abs(BIAS_NO_ADJ-1),
                  DIFF_R2 = R2 - R2_NO_ADJ)
  
  cat(sp_codes[ii], ": ", 
      round(compare_total$BIAS, 2), " (SR); ", 
      round(compare_total$BIAS_NO_ADJ, 2), " (no SR)\n", sep = "")
  
  
}

# Make plots
cowplot::plot_grid(p_fit_sratio_adjustment, p_area_only)
