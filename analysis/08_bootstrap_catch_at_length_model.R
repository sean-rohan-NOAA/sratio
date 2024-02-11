library(sratio)

bootstrap_sample_path <- list.files(here::here("output"), 
                                    recursive = TRUE, 
                                    pattern = "bootstrap_samples_", 
                                    full.names = TRUE)

cal_rmse_df <- read.csv(file = here::here("output", "cal_model_rmse.csv"))

for(ii in 1:length(bootstrap_sample_path)) {
  
  boot_dat <- readRDS(bootstrap_sample_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_sample_path[ii])))
  
  best_family <- cal_rmse_df$name[cal_rmse_df$LOWEST_RMSE & cal_rmse_df$SPECIES_CODE == sp_code]
  
  gam_knots <- cal_rmse_df$gam_knots[cal_rmse_df$LOWEST_RMSE & cal_rmse_df$SPECIES_CODE == sp_code]
  
  bootstrap_df <- sratio::sccal_fit_bootstrap(x = boot_dat,
                                              treatment_order = c(30, 15),
                                              size_col = "SIZE_BIN",
                                              block_col = "MATCHUP",
                                              treatment_col = "TREATMENT",
                                              count_col = "FREQ_EXPANDED",
                                              effort_col = "AREA_SWEPT_KM2",
                                              gam_family = best_family,
                                              k = gam_knots,
                                              n_cores = 4)
  
  # Rename output columns
  bootstrap_df <- bootstrap_df |> 
    select(-AREA_SWEPT_KM2) |>
    dplyr::mutate(SPECIES_CODE = sp_code)
  
  saveRDS(bootstrap_df, file = here::here("output", sp_code,  paste0("sccal_model_bootstrap_results_", sp_code, ".rds")))
  
}


