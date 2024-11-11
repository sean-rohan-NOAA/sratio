# Select best selectivity ratio model

library(sratio)

dat_sratio <- readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds"))

sp_code <- unique(dat_sratio$SPECIES_CODE)

unique_matchups <- unique(dat_sratio$MATCHUP)

pratio_samples <- data.frame()

rmse_results <- data.frame()

for(ii in 1:length(sp_code)) {
  
  spp_lengths <- dplyr::filter(dat_sratio, SPECIES_CODE == sp_code[ii])
  
  # Set knots based on number of length bins, but only use 5 knots for red king crab
  gam_knots <- length(unique(spp_lengths$SIZE_BIN))-4
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(sp_code[ii] %in% c(471, 69322)) {
    gam_knots <- 5
  }
  
  # Run match-up level cross validation
  output_binomial <- 
    sratio_cv(
      model_type = "binomial", 
      count1 = spp_lengths$N_30,
      count2 = spp_lengths$N_15,
      effort1 = spp_lengths$AREA_SWEPT_KM2_30,
      effort2 = spp_lengths$AREA_SWEPT_KM2_15,
      sampling_factor1 = spp_lengths$SAMPLING_FACTOR_30,
      sampling_factor2 = spp_lengths$SAMPLING_FACTOR_15,
      size = spp_lengths$SIZE_BIN,
      block = spp_lengths$MATCHUP,
      k = gam_knots,
      n_cores = 4,
      scale_method = "sv",
      obs_weight_control = 
        list(method = "residuals_by_count", 
             max_count = 200,
             residual_type = "absolute",
             normalize_weight = TRUE)
    )
  
  p_binomial <- output_binomial$cv
  p_binomial$model <- output_binomial$model_settings$model_type
  p_binomial$k <- output_binomial$model_settings$k
  p_binomial$obs_weight_method <- output_binomial$model_settings$obs_weight_control$method
  p_binomial$obs_weight_max_count <- output_binomial$model_settings$obs_weight_control$max_count
  p_binomial$obs_weight_residual_type <- output_binomial$model_settings$obs_weight_control$residual_type
  p_binomial$obs_weight_normalize_weight <- output_binomial$model_settings$obs_weight_control$normalize_weight
  
  output_beta <- 
    sratio_cv(
      model_type = "beta", 
      count1 = spp_lengths$N_30,
      count2 = spp_lengths$N_15,
      effort1 = spp_lengths$AREA_SWEPT_KM2_30,
      effort2 = spp_lengths$AREA_SWEPT_KM2_15,
      sampling_factor1 = spp_lengths$SAMPLING_FACTOR_30,
      sampling_factor2 = spp_lengths$SAMPLING_FACTOR_15,
      size = spp_lengths$SIZE_BIN,
      block = spp_lengths$MATCHUP,
      k = gam_knots,
      n_cores = 4,
      scale_method = "sv",
      obs_weight_control = 
        list(method = "count", 
             max_count = 200,
             residual_type = "none",
             normalize_weight = FALSE)
    )
  
  p_beta <- output_beta$cv
  p_beta$model <- output_beta$model_settings$model_type
  p_beta$k <- output_beta$model_settings$k
  p_beta$obs_weight_method <- output_beta$model_settings$obs_weight_control$method
  p_beta$obs_weight_max_count <- output_beta$model_settings$obs_weight_control$max_count
  p_beta$obs_weight_residual_type <- output_beta$model_settings$obs_weight_control$residual_type
  p_beta$obs_weight_normalize_weight <- output_beta$model_settings$obs_weight_control$normalize_weight
  
  
  # Rename columns to match inputs
  pratio_species <- dplyr::bind_rows(p_binomial, p_beta) |> 
    dplyr::mutate(SPECIES_CODE = sp_code[ii]) |>
    dplyr::select(model,
                  obs_weight_method,
                  obs_weight_max_count,
                  obs_weight_residual_type,
                  obs_weight_normalize_weight,
                  k,
                  SPECIES_CODE,
                  SIZE_BIN = size,
                  MATCHUP = block,
                  N_30 = count1,
                  N_15 = count2,
                  SAMPLING_FACTOR_30 = sampling_factor1,
                  SAMPLING_FACTOR_35 = sampling_factor2,
                  AREA_SWEPT_KM2_30 = effort1,
                  AREA_SWEPT_KM2_15 = effort2,
                  cv_fit,
                  p,
                  s)
                
  pratio_samples <- dplyr::bind_rows(pratio_samples, pratio_species)
  
  # Calculate root mean square error
  rmse_df <- pratio_species |>
    dplyr::group_by(
      SPECIES_CODE, 
      model, 
      k,
      obs_weight_method, 
      obs_weight_max_count, 
      obs_weight_residual_type, 
      obs_weight_normalize_weight
    ) |>
    dplyr::summarise(
      rmse = sqrt(mean((cv_fit - p)^2))
    )
  
  rmse_df$best <- rmse_df$rmse == min(rmse_df$rmse)
  
  rmse_results <- dplyr::bind_rows(rmse_results, rmse_df)
    
}

saveRDS(object = pratio_samples, 
        file = here::here("analysis", "15_30", "output", "pratio_samples.rds"))

write.csv(x = rmse_results, 
          file = here::here("analysis", "15_30", "output", "sratio_model_rmse.csv"), 
          row.names = FALSE)
