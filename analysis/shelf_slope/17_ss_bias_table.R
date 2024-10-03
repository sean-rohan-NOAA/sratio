# Calculate CPUE ratio, bias, mean absolute error, and root mean square error 

library(sratio)

# Load built-in data sets
catch_df <- sratio::data_ss$catch

haul_df <- sratio::data_ss$haul

sp_codes <- sort(unique(catch_df$SPECIES_CODE))

bias_table <- catch_df |>
  dplyr::inner_join(haul_df) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("analysis", "shelf_slope", "output", "n_by_treatment_ss.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_172, AREA_SWEPT_KM2_44) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_172 = WEIGHT_172/AREA_SWEPT_KM2_172,
                CPUE_44 = WEIGHT_44/AREA_SWEPT_KM2_44) |>
  dplyr::mutate(ratio = CPUE_172/CPUE_44,
                log_error = log10(CPUE_172+0.001)-log10(CPUE_44+0.001),
                abs_error = abs(log10(CPUE_172)-log10(CPUE_44)),
                sq_error = (CPUE_172-CPUE_44)^2) |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::summarise(MEAN_RATIO = mean(ratio),
                   BIAS = 10^(mean(log_error)),
                   MAE = 10^mean(abs_error),
                   RMSE = sqrt(mean(sq_error))) |>
  dplyr::mutate(COMMON_NAME = sratio:::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::select(COMMON_NAME, BIAS, MEAN_RATIO, MAE, RMSE)

lines <- c("CPUE comparison between PNE and 83-112 tows\n", "BIAS > 1 = PNE higher\n", "MEAN_RATIO > 1 = PNE higher\n\n\n")
cat(lines, file = here::here("analysis", "shelf_slope", "plots", "bias_table.csv"))

write.table(bias_table, 
            file = here::here("analysis", "shelf_slope", "plots", "bias_table.csv"), 
            append = TRUE, 
            row.names = FALSE, 
            sep = ",")
