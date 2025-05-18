# Select best selectivity ratio model

library(sratio)

dat_sratio <- readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds"))

sp_code <- unique(dat_sratio$SPECIES_CODE)

unique_matchups <- unique(dat_sratio$MATCHUP)

pratio_samples <- data.frame()

rmse_results <- data.frame()

sratio_samples <- data.frame()

binomial_obs_weight_methods <- 
  list(list(method = "count", 
            max_count = 50,
            residual_type = "absolute",
            normalize_weight = FALSE),
       list(method = "residuals_by_count", 
            max_count = 50,
            residual_type = "absolute",
            normalize_weight = FALSE),
       list(method = "none", 
            max_count = Inf,
            residual_type = NA,
            normalize_weight = FALSE)
  )

beta_obs_weight_methods <- 
  list(
    list(method = "none", 
         max_count = 50,
         residual_type = NA,
         normalize_weight = FALSE),
    list(method = "count", 
         max_count = 50,
         residual_type = "none",
         normalize_weight = FALSE)
  )

for(ii in 1:length(sp_code)) {
  
  spp_lengths <- dplyr::filter(dat_sratio, SPECIES_CODE == sp_code[ii]) |>
    dplyr::filter(!(N_15 == 0 & N_30 == 0))
  
  # Set knots based on number of length bins, but only use 5 knots for red king crab
  gam_knots <- length(unique(spp_lengths$SIZE_BIN))-4
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(sp_code[ii] %in% c(471, 69322)) {
    gam_knots <- 5
  }
  
  binomial_fits <- data.frame()
  
  for(jj in 1:length(binomial_obs_weight_methods)) {
    
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
        sratio_type = "absolute",
        obs_weight_control = binomial_obs_weight_methods[[jj]]
      )
    
    p_binomial <- output_binomial$cv
    p_binomial$model <- output_binomial$model_settings$model_type
    p_binomial$k <- output_binomial$model_settings$k
    p_binomial$obs_weight_method <- output_binomial$model_settings$obs_weight_control$method
    p_binomial$obs_weight_max_count <- output_binomial$model_settings$obs_weight_control$max_count
    p_binomial$obs_weight_residual_type <- output_binomial$model_settings$obs_weight_control$residual_type
    p_binomial$obs_weight_normalize_weight <- output_binomial$model_settings$obs_weight_control$normalize_weight
    
    binomial_fits <- dplyr::bind_rows(binomial_fits, p_binomial)
    
    # Add calculated selectivity ratios to a data.frame to be used as a built-in data set.
    if(jj == 1) {
      sratio_samples <- output_binomial$data |>
        dplyr::mutate(species_code = sp_code[ii]) |>
        dplyr::bind_rows(sratio_samples)
    }
    
    
  }
  
  beta_fits <- data.frame()
  
  for(kk in 1:length(beta_obs_weight_methods)) {
    
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
        sratio_type = "absolute",
        obs_weight_control = beta_obs_weight_methods[[kk]]
      )
    
    p_beta <- output_beta$cv
    p_beta$model <- output_beta$model_settings$model_type
    p_beta$k <- output_beta$model_settings$k
    p_beta$obs_weight_method <- output_beta$model_settings$obs_weight_control$method
    p_beta$obs_weight_max_count <- output_beta$model_settings$obs_weight_control$max_count
    p_beta$obs_weight_residual_type <- output_beta$model_settings$obs_weight_control$residual_type
    p_beta$obs_weight_normalize_weight <- output_beta$model_settings$obs_weight_control$normalize_weight
    
    beta_fits <- dplyr::bind_rows(beta_fits, p_beta)
    
  }
  
  # Rename columns to match inputs
  pratio_species <- 
    dplyr::bind_rows(binomial_fits, beta_fits) |> 
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
                  SAMPLING_FACTOR_15 = sampling_factor2,
                  AREA_SWEPT_KM2_30 = effort1,
                  AREA_SWEPT_KM2_15 = effort2,
                  p,
                  s,
                  p_fit,
                  s_fit)
  
  # Calculate root mean square error for proportions
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
      rmse = sqrt(mean((p_fit - p)^2))
    )

  rmse_df$best <- rmse_df$rmse == min(rmse_df$rmse)
  
  rmse_results <- dplyr::bind_rows(rmse_results, rmse_df)
  
  pratio_species$FIT_N_30 <- 
    sratio_predict_count(
      est_count = "count1",
      count1 = pratio_species$N_30, 
      count2 = pratio_species$N_15, 
      sampling_factor1 = pratio_species$SAMPLING_FACTOR_30, 
      sampling_factor2 = pratio_species$SAMPLING_FACTOR_15, 
      effort1 = pratio_species$AREA_SWEPT_KM2_30, 
      effort2 = pratio_species$AREA_SWEPT_KM2_15, 
      s12 = pratio_species$s_fit
    )
  
  pratio_species$FIT_N_15 <- 
    sratio_predict_count(
      est_count = "count2",
      count1 = pratio_species$N_30, 
      count2 = pratio_species$N_15, 
      sampling_factor1 = pratio_species$SAMPLING_FACTOR_30, 
      sampling_factor2 = pratio_species$SAMPLING_FACTOR_15, 
      effort1 = pratio_species$AREA_SWEPT_KM2_30, 
      effort2 = pratio_species$AREA_SWEPT_KM2_15, 
      s12 = pratio_species$s_fit
    )
  
  pratio_samples <- dplyr::bind_rows(pratio_samples, pratio_species)
    
}

save(sratio_samples, file = here::here("analysis", "15_30", "data", "sratio_samples_1530.rda"))

saveRDS(object = pratio_samples, 
        file = here::here("analysis", "15_30", "output", "pratio_samples.rds"))

write.csv(x = rmse_results, 
          file = here::here("analysis", "15_30", "output", "sratio_model_rmse.csv"), 
          row.names = FALSE)
