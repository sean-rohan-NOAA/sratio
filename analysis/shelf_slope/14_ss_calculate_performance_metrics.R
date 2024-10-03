# Compare predicted catch in 30 minute hauls to observed catch-at-length
# Based on selectivity ratio and selectivity condition on catch-at-length models.

library(sratio)

sp_codes <- c(21740, 21720, 10110, 10112, 10115, 10130, 471,  68580)

# Use the median or mean of bootstrap samples?
central_fun <- mean

# Tendency measure - Calculate the median or mean of performance measures?
calc_type <- mean

model_method <- c("sratio", "sccal")

for(jj in 1:length(model_method)) {
  
  # Minimum catch-at-length for 30 minute tows to be used for performance metrics
  min_n <- 5
  
  for(ii in 1:length(sp_codes)) {
    
    # Load catch-at-length and effort data
    size_172 <- readRDS(here::here("analysis", "shelf_slope", "output", "catch_at_length_1530.rds")) |>
      dplyr::filter(TREATMENT == 172, SPECIES_CODE == sp_codes[ii]) |>
      dplyr::ungroup()
    
    size_44 <- readRDS(here::here("analysis", "shelf_slope", "output", "catch_at_length_1530.rds")) |>
      dplyr::filter(TREATMENT == 44, SPECIES_CODE == sp_codes[ii]) |>
      dplyr::ungroup()
    
    # Load bootstrap results from best-fitting model
    bootstrap_results_path <- list.files(here::here("analysis", "shelf_slope", "output"), 
                                         recursive = TRUE, 
                                         pattern = paste0(model_method[jj], "_bootstrap_results_"), 
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
    
    size_44 <- dplyr::select(size_172, MATCHUP, AREA_SWEPT_KM_30 = AREA_SWEPT_KM2) |> 
      unique() |>
      dplyr::inner_join(size_44, by = join_by(MATCHUP)) |>
      dplyr::inner_join(bootstrap_results, by = join_by(SIZE_BIN))
    
    # Predict mean catch-at-length, adjusting for area swept and selectivity
    area_swept_ratio <- size_44$AREA_SWEPT_KM_30/size_44$AREA_SWEPT_KM2
    
    size_44$PREDICTED_FREQUENCY <- size_44$FREQ_EXPANDED * area_swept_ratio * size_44$s12
    
    # Predict mean catch-at-length adjusting for area swept but not selectivity
    size_44$PREDICTED_FREQUENCY_NO_ADJ <- size_44$FREQ_EXPANDED * area_swept_ratio
    
    compare_fit <- size_44 |> 
      dplyr::select(MATCHUP, 
                    PREDICTED_FREQUENCY, 
                    PREDICTED_FREQUENCY_NO_ADJ, 
                    SIZE_BIN, 
                    SPECIES_CODE) |>
      dplyr::full_join(size_172, 
                       by = join_by(MATCHUP, SIZE_BIN, SPECIES_CODE)) |>
      dplyr::mutate(PREDICTED_FREQUENCY = if_else(is.na(PREDICTED_FREQUENCY), 0, PREDICTED_FREQUENCY),
                    PREDICTED_FREQUENCY_NO_ADJ = if_else(is.na(PREDICTED_FREQUENCY_NO_ADJ), 0, PREDICTED_FREQUENCY_NO_ADJ),
                    FREQ_EXPANDED = if_else(is.na(FREQ_EXPANDED), 0, FREQ_EXPANDED)) |>
      dplyr::filter(FREQ_EXPANDED >= min_n) # Filter records with insufficient samples
    
    # Calculate performance metrics
    compare_by_size <- compare_fit |> 
      dplyr::group_by(SIZE_BIN, SPECIES_CODE) |>
      dplyr::summarise(N_HAUL = n(),
                       MAE = calc_mae(est = PREDICTED_FREQUENCY, 
                                      obs = FREQ_EXPANDED, 
                                      const = 1e-3,
                                      type = median),
                       MRE = suppressWarnings(calc_mre(est = PREDICTED_FREQUENCY, 
                                                       obs = FREQ_EXPANDED,
                                                       type = calc_type)),
                       BIAS = calc_bias(est = PREDICTED_FREQUENCY, 
                                        obs = FREQ_EXPANDED, 
                                        const = 1e-3, 
                                        type = calc_type),
                       R2 = suppressWarnings(cor(PREDICTED_FREQUENCY, 
                                                 FREQ_EXPANDED, 
                                                 use = "complete.obs")^2),
                       RMSE = calc_rmse(est = PREDICTED_FREQUENCY, 
                                        obs = FREQ_EXPANDED,
                                        type = calc_type),
                       MAE_NO_ADJ = calc_mae(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = FREQ_EXPANDED, 
                                             const = 1e-3,
                                             type = calc_type),
                       MRE_NO_ADJ = suppressWarnings(calc_mre(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                                              obs = FREQ_EXPANDED,
                                                              type = calc_type)),
                       BIAS_NO_ADJ = calc_bias(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                               obs = FREQ_EXPANDED, 
                                               const = 1e-3, 
                                               type = calc_type),
                       R2_NO_ADJ = suppressWarnings(cor(PREDICTED_FREQUENCY_NO_ADJ, 
                                                        FREQ_EXPANDED, 
                                                        use = "complete.obs")^2),
                       RMSE_NO_ADJ = calc_rmse(est = PREDICTED_FREQUENCY_NO_ADJ, 
                                               obs = FREQ_EXPANDED,
                                               type = calc_type),
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
                                      const = 1e-3,
                                      type = calc_type),
                       MRE = calc_mre(est = TOTAL_PREDICTED_FREQUENCY, 
                                      obs = TOTAL_FREQ_EXPANDED,
                                      type = calc_type),
                       BIAS = calc_bias(est = TOTAL_PREDICTED_FREQUENCY, 
                                        obs = TOTAL_FREQ_EXPANDED, 
                                        const = 1e-3,
                                        type = calc_type),
                       R2 = cor(TOTAL_PREDICTED_FREQUENCY, 
                                TOTAL_FREQ_EXPANDED, 
                                use = "complete.obs")^2,
                       RMSE = calc_rmse(est = TOTAL_PREDICTED_FREQUENCY, 
                                        obs = TOTAL_FREQ_EXPANDED,
                                        type = calc_type),
                       MAE_NO_ADJ = calc_mae(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = TOTAL_FREQ_EXPANDED, 
                                             const = 1e-3,
                                             type = calc_type),
                       MRE_NO_ADJ = calc_mre(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                             obs = TOTAL_FREQ_EXPANDED),
                       BIAS_NO_ADJ = calc_bias(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                               obs = TOTAL_FREQ_EXPANDED, 
                                               const = 1e-3,
                                               type = calc_type),
                       R2_NO_ADJ = cor(TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                       TOTAL_FREQ_EXPANDED, 
                                       use = "complete.obs")^2,
                       RMSE_NO_ADJ = calc_rmse(est = TOTAL_PREDICTED_FREQUENCY_NO_ADJ, 
                                               obs = TOTAL_FREQ_EXPANDED,
                                               type = calc_type),
                       .groups = "keep") |>
      dplyr::mutate(DIFF_RMSE = RMSE - RMSE_NO_ADJ,
                    DIFF_MAE = MAE - MAE_NO_ADJ,
                    DIFF_MRE = abs(MRE) - abs(MRE_NO_ADJ),
                    DIFF_BIAS = abs(BIAS-1) - abs(BIAS_NO_ADJ-1),
                    DIFF_R2 = R2 - R2_NO_ADJ)
    
    saveRDS(object = list(species_code = sp_codes[ii],
                          compare_fit = compare_fit,
                          performance_metrics_by_size = compare_by_size,
                          performance_metrics = compare_total,
                          performance_metrics_by_haul = compare_total_by_haul),
            file = here::here("analysis", "shelf_slope", "output", paste0("model_performance_", model_method[jj], "_", sp_codes[ii], ".rds")))
    
    cat(sp_codes[ii], ": ", 
        round(compare_total$BIAS, 2), " (SR); ", 
        round(compare_total$BIAS_NO_ADJ, 2), " (no SR)\n", sep = "")
    
  }
  
}
