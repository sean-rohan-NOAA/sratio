library(sratio)

bootstrap_sample_path <- list.files(here::here("analysis", "shelf_slope", "output"), 
                                    recursive = TRUE, 
                                    pattern = "bootstrap_samples_", 
                                    full.names = TRUE)

sratio_rmse_df <- read.csv(file = here::here("analysis", "shelf_slope", "output", "sratio_model_rmse.csv"))

for(ii in 1:length(bootstrap_sample_path)) {
  
  boot_dat <- readRDS(bootstrap_sample_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_sample_path[ii])))
  
  # Set model type (binomial or beta)
  best_model <- sratio_rmse_df$model[sratio_rmse_df$best & sratio_rmse_df$SPECIES_CODE == sp_code]
  
  gam_knots <- sratio_rmse_df$gam_knots[sratio_rmse_df$best & sratio_rmse_df$SPECIES_CODE == sp_code]
  
  # Run two-stage bootstrap (stage 1: by matchup, stage 2: by sample within a matchup)
  bootstrap_df <- sratio::sratio_fit_bootstrap(x = boot_dat,
                                               treatment_order = c(172, 44),
                                               size_col = "SIZE_BIN",
                                               block_col = "MATCHUP",
                                               treatment_col = "TREATMENT",
                                               count_col = "FREQ_EXPANDED",
                                               effort_col = "AREA_SWEPT_KM2",
                                               gam_family = best_model,
                                               k = gam_knots,
                                               scale_method = "sv",
                                               n_cores = 4)
  bootstrap_df$SPECIES_CODE <- sp_code
  
  saveRDS(bootstrap_df, 
          here::here("analysis", "shelf_slope", "output", 
                     sp_code, 
                     paste0("sratio_bootstrap_results_", sp_code, ".rds")))
  
}


bootstrap_results_path <- list.files(here::here("analysis", "shelf_slope", "output"), 
                                     recursive = TRUE, 
                                     pattern = "sratio_bootstrap_results_", 
                                     full.names = TRUE)

bootstrap_samples <- data.frame()

for(kk in 1:length(bootstrap_results_path)) {
  
  bootstrap_samples <- dplyr::bind_rows(bootstrap_samples,
                                        readRDS(bootstrap_results_path[kk]))
  
}

bootstrap_samples <- bootstrap_samples |>
  dplyr::select(SPECIES_CODE, SIZE_BIN, p12, s12)


saveRDS(object = bootstrap_samples,
        file = here::here("analysis", "shelf_slope", "output", "shelf_slope_sratio_bootstrap.rds"))

