# Calculate CPUE ratio, bias, mean absolute error, and root mean square error 

library(sratio)

seed <- 1729730909

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
  dplyr::reframe(BIAS = bootstrap_mean_bias(CPUE_30, CPUE_15, n_samples = 10000, add_constant = 1, scale = "log10", seed = seed)) |>
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

saveRDS(bias_quantiles, 
        file = here::here("analysis", "15_30", "output", "cpue_bias_bootstrap_quantiles.rds"))


lines <- c("CPUE comparison between 15 and 30 minute tows\n", "BIAS > 1 = 30 minutes higher\n\n\n")
cat(lines, file = here::here("analysis", "15_30", "plots", "total_cpue_fit", "bias_table.csv"))

write.table(bias_table, 
            file = here::here("analysis", "15_30", "plots", "total_cpue_fit", "bias_table.csv"), 
            append = TRUE, 
            row.names = FALSE, 
            sep = ",")
