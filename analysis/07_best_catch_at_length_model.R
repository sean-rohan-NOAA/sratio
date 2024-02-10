# Select best catch-at-length model

library(sratio)

dat_binned <- readRDS(file = here::here("output", "catch_at_length_1530.rds"))
# dat_binned$FREQ_EXPANDED <- round(dat_binned$FREQ_EXPANDED)

rmse_df <- data.frame()

species_codes <- unique(dat_binned$SPECIES_CODE)

for(ii in 1:length(species_codes)) {
  
  spp_lengths <- dplyr::filter(dat_binned, SPECIES_CODE == species_codes[ii])
  spp_lengths$dummy_var <- 1
  
  gam_knots <- (length(unique(spp_lengths $SIZE_BIN))-1)-3
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(species_codes[ii] == 69322) {
    gam_knots <- 5
  }
  
  # Setup four clusters and folds for each matchups
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  folds <- caret::groupKFold(group = interaction(spp_lengths$MATCHUP, spp_lengths$TREATMENT))
  
  cv_results <- foreach::foreach(fold = folds, .packages = "mgcv") %dopar% {
    
    training_df <- spp_lengths[fold, ]
    validation_df <- spp_lengths[-fold, ]
    validation_df$dummy_var <- 0
    
    # Add in dummy station variable for predictions, to be added back in for output
    out_matchup <- validation_df$MATCHUP[1]
    validation_df$MATCHUP <- training_df$MATCHUP[1]
    
    mod_tw <- mgcv::gam(formula = FREQ_EXPANDED ~ s(SIZE_BIN, k = gam_knots, bs = "tp", by = TREATMENT) + 
                          s(MATCHUP, 
                            by = dummy_var,
                            bs = "re") + 
                          offset(I(log(AREA_SWEPT_KM2))), 
                        data = training_df,
                        family = tw(link = "log"))
    
    mod_nb <- mgcv::gam(formula = FREQ_EXPANDED ~ s(SIZE_BIN, k = gam_knots, bs = "tp", by = TREATMENT) + 
                          s(MATCHUP, 
                            by = dummy_var,
                            bs = "re") + 
                          offset(I(log(AREA_SWEPT_KM2))), 
                        data = training_df,
                        family = nb(link = "log"))
    
    mod_poisson <- mgcv::gam(formula = FREQ_EXPANDED ~ s(SIZE_BIN, k = gam_knots, bs = "tp", by = TREATMENT) + 
                               s(MATCHUP, 
                                 by = dummy_var,
                                 bs = "re") + 
                               offset(I(log(AREA_SWEPT_KM2))), 
                             data = training_df,
                             family = poisson(link = "log"))
    
    validation_df$tw <- predict(mod_tw, newdata = validation_df, type = "response")
    validation_df$nb <- predict(mod_nb, newdata = validation_df, type = "response")
    validation_df$poisson <- predict(mod_poisson, newdata = validation_df, type = "response")
    
    # Reset matchup and dummy variable for fitting final models
    validation_df$MATCHUP <- out_matchup
    validation_df$dummy_var <- 1
    
    return(validation_df)
  }
  
  results <- do.call("rbind", cv_results)
  
  doParallel::stopImplicitCluster()
  
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

write.csv(rmse_df, here::here("output", "cal_model_rmse.csv"), row.names = FALSE)
