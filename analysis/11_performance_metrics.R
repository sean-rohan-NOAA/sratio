# Calculate CPUE ratio, bias, mean absolute error, and root mean square error

library(sratio)

# Load built-in data sets
catch_df <- sratio::data_1530$catch |>
  dplyr::filter(CRUISE %in% use_cruises)

haul_df <- sratio::data_1530$haul |>
  dplyr::filter(CRUISE %in% use_cruises)

sp_codes <- sort(unique(catch_df$SPECIES_CODE))

bias_table <- catch_df |>
  dplyr::inner_join(haul_df) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("output", "n_by_treatment_1530.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
                CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
  dplyr::mutate(ratio = CPUE_15/CPUE_30,
                log_error = log10(CPUE_15+0.001)-log10(CPUE_30+0.001),
                abs_error = abs(log10(CPUE_15)-log10(CPUE_30)),
                sq_error = (CPUE_15-CPUE_30)^2) |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::summarise(MEAN_RATIO = mean(ratio),
                   BIAS = 10^(mean(log_error)),
                   MAE = 10^mean(abs_error),
                   RMSE = sqrt(mean(sq_error))) |>
  dplyr::mutate(COMMON_NAME = sratio:::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::select(COMMON_NAME, BIAS, MEAN_RATIO, MAE, RMSE)

lines <- c("CPUE comparison between 15 and 30 minute tows\n", "BIAS > 1 = 15 minutes higher\n", "MEAN_RATIO > 1 = 15 minutes higher\n\n\n")
cat(lines, file = here::here("plots", "bias_table.csv"))

write.table(bias_table, 
            file = here::here("plots", "bias_table.csv"), 
            append = TRUE, 
            row.names = FALSE, 
            sep = ",")
