library(sratio)

bootstrap_sample_path <- list.files(here::here("analysis", "15_30", "output"), 
                                    recursive = TRUE, 
                                    pattern = "bootstrap_samples_", 
                                    full.names = TRUE)

discrete_rmse_df <- read.csv(file = here::here("analysis", "15_30", "output", "discrete_model_rmse.csv"))

for(ii in 1:length(bootstrap_sample_path)) {
  
  boot_dat <- readRDS(bootstrap_sample_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_sample_path[ii])))
  
  # Set model type (binomial or beta)
  best_model <- discrete_rmse_df$model[discrete_rmse_df$best & discrete_rmse_df$SPECIES_CODE == sp_code]
  
  gam_knots <- discrete_rmse_df$gam_knots[discrete_rmse_df$best & discrete_rmse_df$SPECIES_CODE == sp_code]
  
  # Run two-stage bootstrap (stage 1: by matchup, stage 2: by sample within a matchup)
  bootstrap_df <- sratio::discrete_fit_bootstrap(x = boot_dat,
                                                 treatment_order = c(30, 15),
                                                 size_col = "SIZE_BIN",
                                                 block_col = "MATCHUP",
                                                 treatment_col = "TREATMENT",
                                                 count_col = "FREQ_EXPANDED",
                                                 effort_col = "AREA_SWEPT_KM2",
                                                 gam_family = best_model,
                                                 scale_method = "sv",
                                                 n_cores = 4)
  bootstrap_df$SPECIES_CODE <- sp_code
  
  saveRDS(bootstrap_df, here::here("analysis", "15_30", "output", 
                                   sp_code, 
                                   paste0("discrete_bootstrap_results_", sp_code, ".rds")))
  
}
