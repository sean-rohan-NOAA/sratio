# Select best catch-at-length model

library(sratio)

dat_binned <- readRDS(file = here::here("analysis", "15_30", "output", "catch_at_length_1530.rds"))
# dat_binned$FREQ_EXPANDED <- round(dat_binned$FREQ_EXPANDED)

rmse_df <- data.frame()

species_codes <- unique(dat_binned$SPECIES_CODE)

for(ii in 1:length(species_codes)) {
  
  spp_lengths <- dplyr::filter(dat_binned, SPECIES_CODE == species_codes[ii]) |>
    dplyr::mutate(MATCHUP = factor(MATCHUP))
  spp_lengths$dummy_var <- 1
  
  # Set GAM knots
  gam_knots <- (length(unique(spp_lengths $SIZE_BIN))-1)-3
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(species_codes[ii] == 69322) {
    gam_knots <- 5
  }
  
  # Run haul level cross-validation on selectivity conditioned on catch-at-length models
  results <- sratio::sccal_cv(count1 = spp_lengths$FREQ_EXPANDED[spp_lengths$TREATMENT == 30], 
                              count2 = spp_lengths$FREQ_EXPANDED[spp_lengths$TREATMENT == 15], 
                              effort1 = spp_lengths$AREA_SWEPT_KM2[spp_lengths$TREATMENT == 30], 
                              effort2 = spp_lengths$AREA_SWEPT_KM2[spp_lengths$TREATMENT == 15], 
                              size1 = spp_lengths$SIZE_BIN[spp_lengths$TREATMENT == 30], 
                              size2 = spp_lengths$SIZE_BIN[spp_lengths$TREATMENT == 15], 
                              block1 = spp_lengths$MATCHUP[spp_lengths$TREATMENT == 30], 
                              block2 = spp_lengths$MATCHUP[spp_lengths$TREATMENT == 15],
                              treatment_name1 = 30,
                              treatment_name2 = 15,
                              k = gam_knots, 
                              n_cores = 4)
  
  # Rename columns
  results <- results |> 
    dplyr::rename(FREQ_EXPANDED = count,
                  SIZE_BIN = size,
                  AREA_SWEPT_KM2 = effort,
                  MATCHUP = block,
                  TREATMENT = treatment)
  
  rmse_df <- results |>
    tidyr::pivot_longer(cols = c("tw", "nb", "poisson")) |>
    dplyr::mutate(resid = (value-FREQ_EXPANDED)^2) |>
    dplyr::group_by(name) |>
    dplyr::summarize(rmse = sqrt(mean(resid, na.rm = TRUE))) |>
    dplyr::mutate(SPECIES_CODE = species_codes[ii],
                  LOWEST_RMSE = min(rmse) == rmse,
                  DELTA_RMSE = rmse-min(rmse), 
                  gam_knots = gam_knots) |>
    as.data.frame() |>
    dplyr::bind_rows(rmse_df)
  
}

write.csv(rmse_df, here::here("analysis", "15_30", "output", "cal_model_rmse.csv"), 
          row.names = FALSE)
