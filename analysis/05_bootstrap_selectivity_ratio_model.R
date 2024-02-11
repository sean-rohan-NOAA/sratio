library(sratio)

n_cores <- 4

# Function to format bootstrap data frame
make_sratio_df <- function(x) {
  
  lengths <- x |>
    dplyr::mutate(TREATMENT_COL = paste0("N_", TREATMENT)) |>
    dplyr::select(MATCHUP, SIZE_BIN, FREQ_EXPANDED, TREATMENT_COL) |>
    tidyr::pivot_wider(names_from = "TREATMENT_COL", values_from = "FREQ_EXPANDED", values_fill = 0) |>
    as.data.frame()
  
  effort <- x |>
    dplyr::mutate(AREA_SWEPT_COL = paste0("AREA_SWEPT_KM2_", TREATMENT)) |>
    dplyr::select(AREA_SWEPT_COL, AREA_SWEPT_KM2, MATCHUP) |>
    unique() |>
    tidyr::pivot_wider(names_from = "AREA_SWEPT_COL", values_from = "AREA_SWEPT_KM2", values_fill = 0) |>
    as.data.frame()
  
  combined <- dplyr::inner_join(lengths, effort, by = "MATCHUP")
  
  combined$p <- suppressMessages(sratio::selectivity_ratio(count1 = combined$N_30, 
                                                           count2 = combined$N_15, 
                                                           effort1 = combined$AREA_SWEPT_KM2_30, 
                                                           effort2 = combined$AREA_SWEPT_KM2_15)$p12)
  
  combined$p_scaled <- sratio::scale_for_betareg(combined$p, method = "sv")
  combined$dummy_var <- 1
  
  return(combined)
  
}

bootstrap_sample_path <- list.files(here::here("output"), 
                                    recursive = TRUE, 
                                    pattern = "bootstrap_samples_", 
                                    full.names = TRUE)

sratio_rmse_df <- read.csv(file = here::here("output", "sratio_model_rmse.csv"))

for(ii in 1:length(bootstrap_sample_path)) {
  
  boot_dat <- readRDS(bootstrap_sample_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_sample_path[ii])))
  
  best_model <- sratio_rmse_df$model[sratio_rmse_df$best & sratio_rmse_df$SPECIES_CODE == sp_code]
  
  gam_knots <- sratio_rmse_df$gam_knots[sratio_rmse_df$best & sratio_rmse_df$SPECIES_CODE == sp_code]
  
  test <- sratio_bootstrap(x = boot_dat, 
                   size_col = "SIZE_BIN", 
                   block_col = "MATCHUP", 
                   treatment_col = "TREATMENT", 
                   count_col = "FREQ_EXPANDED", 
                   effort_col = "AREA_SWEPT_KM2", 
                   gam_family = best_model, 
                   k = gam_knots, 
                   scale_method = "sv", 
                   n_cores = 4)
  
  # gam_family <- switch(best_model,
  #                      "Binomial" = binomial(link = "logit"),
  #                      "Beta" = betar(link = "logit"))
  # 
  # boot_dat <- lapply(boot_dat, make_sratio_df)
  # 
  # lengths <- seq(min(unlist(lapply(boot_dat, FUN = function(x) {min(x$SIZE_BIN)}))), 
  #                max(unlist(lapply(boot_dat, FUN = function(x) {max(x$SIZE_BIN)}))),
  #                by = 1)
  
  # cl <- parallel::makeCluster(n_cores)
  # doParallel::registerDoParallel(cl)
  # 
  # bootstrap_output <- foreach::foreach(iter = 1:length(boot_dat), .packages = c("mgcv", "dplyr")) %dopar% {
  #   
  #   boot_df <- boot_dat[[iter]] |>
  #     dplyr::mutate(MATCHUP = factor(MATCHUP))
  #   
  #   gam_formula <- switch(best_model,
  #                         "Binomial" =  p ~ s(SIZE_BIN, k = gam_knots, bs = 'tp') + s(MATCHUP, bs = 're', by = dummy_var),
  #                         "Beta" =  p_scaled ~ s(SIZE_BIN, k = gam_knots, bs = 'tp') + s(MATCHUP, bs = 're', by = dummy_var))
  #   
  #   model <- mgcv::gam(formula = gam_formula,
  #                      family = gam_family,
  #                      data = boot_df)
  #   
  #   fit_df <- data.frame(SIZE_BIN = lengths,
  #                        MATCHUP = boot_dat[[iter]]$MATCHUP[1],
  #                        dummy_var = 0) # random effects off
  #   
  #   fit_df$fit <- predict(model, newdata = fit_df, type = "link")
  #   
  #   return(fit_df)
  #   
  # }
  # 
  # doParallel::stopImplicitCluster()
  # 
  # bootstrap_df <- do.call("rbind", bootstrap_output)
  
  saveRDS(bootstrap_df, here::here("output", sp_code, paste0("sratio_bootstrap_results_", sp_code, ".rds")))
  
}
