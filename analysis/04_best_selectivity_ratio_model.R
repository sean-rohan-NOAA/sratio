# Select best selectivity ratio model

library(sratio)

dat_sratio <- readRDS(file = here::here("output", "n_by_treatment_1530.rds"))

species_codes <- unique(dat_sratio$SPECIES_CODE)

#temp species drop
species_codes <- species_codes[-which(species_codes == 68580)]

unique_matchups <- unique(dat_sratio$MATCHUP)

n_cores <- 4

rmse_df <- data.frame()

for(ii in 1:length(species_codes)) {
  
  pratio_df <- data.frame()
  
  spp_lengths <- dplyr::filter(dat_sratio, SPECIES_CODE == species_codes[ii])
  
  for(jj in 1:length(unique_matchups)) {
    
    sel_length <- dplyr::filter(spp_lengths, MATCHUP == unique_matchups[jj])
    
    # Calculate absolute selectivity ratio
    sel_length$p <- suppressMessages(
      selectivity_ratio(count1 = sel_length$N_30, 
                        count2 = sel_length$N_15, 
                        effort1 = sel_length$AREA_SWEPT_KM2_30, 
                        effort2 = sel_length$AREA_SWEPT_KM2_15)$p12
    )
    
    pratio_df <- dplyr::bind_rows(pratio_df, sel_length)
    
  }
  
  pratio_df <- dplyr::filter(pratio_df, !is.na(p))
  pratio_df$p_scaled <- sratio::scale_for_betareg(pratio_df$p, method = "sv")
  pratio_df$dummy_var <- 1
  
  # Set knots based on number of length bins, but only use 5 knots for red king crab
  gam_knots <- (length(unique(pratio_df$SIZE_BIN))-1)-3
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(species_codes[ii] == 69322) {
    gam_knots <- 5
  }
  
  # Setup four clusters and folds for each matchups
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  folds <- caret::groupKFold(group = pratio_df$MATCHUP)
  
  cv_results <- foreach::foreach(fold = folds, .packages = c("mgcv", "dplyr")) %dopar% {
    
    training_df <- pratio_df[fold, ]
    validation_df <- pratio_df[-fold, ]
    validation_df$dummy_var <- 0
    
    # Add in dummy station variable for predictions, to be added back in for output
    out_matchup <- validation_df$MATCHUP[1]
    validation_df$MATCHUP <- training_df$MATCHUP[1]
    
    gam_logit <- mgcv::gam(p_scaled ~ s(SIZE_BIN, bs = "tp", k = gam_knots) + s(MATCHUP, bs = "re", by = dummy_var),
                           data = training_df |>
                             dplyr::mutate(MATCHUP = factor(MATCHUP)),
                           family = binomial(link = "logit"))
    
    gam_beta <- mgcv::gam(p_scaled ~ s(SIZE_BIN, bs = "tp", k = gam_knots) + s(MATCHUP, bs = "re", by = dummy_var),
                          data = training_df |>
                            dplyr::mutate(MATCHUP = factor(MATCHUP)),
                          family = mgcv::betar(link = "logit"))
    
    fitted_logit <- predict(gam_logit, newdata = validation_df, type = "response")
    fitted_beta <- predict(gam_beta, newdata = validation_df, type = "response")
    
    validation_df$cv_fit_logit <- fitted_logit
    validation_df$cv_fit_beta <- fitted_beta
    
    # Reset matchup and dummy variable for fitting final models
    validation_df$MATCHUP <- out_matchup
    validation_df$dummy_var <- 1
    
    return(validation_df)
  }
  
  doParallel::stopImplicitCluster()
  
  pratio_df <- do.call("rbind", cv_results)

  rmse_df <- rmse_df |>
    dplyr::bind_rows(
      data.frame(model = c("Binomial", "Beta"), 
                 rmse = c(sqrt(mean((pratio_df$cv_fit_logit-pratio_df$p_scaled)^2)),
                          sqrt(mean((pratio_df$cv_fit_beta-pratio_df$p_scaled)^2))),
                 SPECIES_CODE = species_codes[ii],
                 gam_knots = gam_knots) |>
        dplyr::mutate(best = rmse == min(rmse))
    )
}

write.csv(rmse_df, here::here("output", "sratio_model_rmse.csv"), row.names = FALSE)
