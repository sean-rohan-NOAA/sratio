binom_cv <- function(count1, count2, effort1, effort2, size, block, k = NULL, scale_method = "sv", n_cores = 1) {
  
  ####
  dat_sratio <- readRDS(file = here::here("output", "n_by_treatment_1530.rds"))
  
  sp_code <- unique(dat_sratio$SPECIES_CODE)
  
  #temp species drop
  if(!any(use_cruises %in% c(199501, 199801))) {
    sp_code <- sp_code[-which(sp_code == 68580)]
  }
  
  unique_matchups <- unique(dat_sratio$MATCHUP)
  
  # n_cores <- 4
  
  rmse_df <- data.frame()
  
  pratio_samples <- data.frame()
  
  ii <- 2
    
    pratio_df <- data.frame()
    
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
    count1 = spp_lengths$N_30
    count2 = spp_lengths$N_15
    effort1 = spp_lengths$AREA_SWEPT_KM2_30
    effort2 = spp_lengths$AREA_SWEPT_KM2_15
    size = spp_lengths$SIZE_BIN
    block = spp_lengths$MATCHUP
    k = gam_knots
    n_cores = 4
    scale_method = "sv"
    
    ###########
  
  unique_blocks <- unique(block)
  
  p <- numeric(length = length(count1))
  s <- numeric(length = length(count1))
  
  for(jj in 1:length(unique_blocks)) {
    
    s_ratio_df <- suppressMessages(
      selectivity_ratio(count1 = count1[block == unique_blocks[jj]], 
                        count2 = count2[block == unique_blocks[jj]],
                        effort1 = effort1[block == unique_blocks[jj]], 
                        effort2 = effort2[block == unique_blocks[jj]] 
      )
    )
    
    p[block == unique_blocks[jj]] <- s_ratio_df$p12
    s[block == unique_blocks[jj]] <- s_ratio_df$s12
    
  }
  
  # Remove NA catch comparison rate
  count1 <- count1[!is.na(p)]
  count2 <- count2[!is.na(p)]
  effort1 <- effort1[!is.na(p)]
  effort2 <- effort2[!is.na(p)]
  size <- size[!is.na(p)]
  block <- block[!is.na(p)]
  s <- p[!is.na(p)]
  p <- p[!is.na(p)]
  
  # Scale p to fall within (the interval supported by (0,1] for beta regression
  p_scaled <- scale_for_betareg(p, method = scale_method)
  
  # Set k if not provided
  if(is.null(k)) {
    k <- min(c(8, (length(unique(size))-3)))
  }
  
  model_df <- data.frame(block = factor(block),
                         size = size,
                         count1 = count1,
                         count2 = count2,
                         effort1 = effort1,
                         effort2 = effort2,
                         p = p,
                         p_scaled = p_scaled,
                         s = s,
                         dummy_var = 1)
  
  # Setup four clusters and folds for each block
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  folds <- caret::groupKFold(group = block)
  
  ###
  fold <- folds[[1]]
  ###
  
  cv_results <- foreach::foreach(fold = folds, .packages = c("mgcv", "dplyr")) %dopar% {
    
    training_df <- model_df[fold, ]
    validation_df <- model_df[-fold, ]
    validation_df$dummy_var <- 0
    
    # Add in dummy station variable for predictions, to be added back in for output
    out_matchup <- validation_df$block[1]
    validation_df$block <- training_df$block[1]
    
    gam_logit <- mgcv::gam(cbind(count1, count2) 
                           ~ s(size, bs = "tp", k = k) + s(block, bs = "re", by = dummy_var) + 
                             offset(I(log(effort1/effort2))),
                           data = training_df,
                           family = binomial())
    
    gam_beta <- mgcv::gam(cbind(count1, count2) 
                          ~ s(size, bs = "tp", k = k) + s(block, bs = "re", by = dummy_var) + 
                            offset(I(log(effort1/effort2))),
                          data = training_df,
                          family = mgcv::betar())
    
    fitted_logit <- predict(gam_logit, newdata = validation_df, type = "response")
    fitted_beta <- predict(gam_beta, newdata = validation_df, type = "response")
    
    validation_df$cv_fit_logit <- fitted_logit
    validation_df$cv_fit_beta <- fitted_beta
    
    # Reset matchup and dummy variable for fitting final models
    validation_df$block <- out_matchup
    validation_df$dummy_var <- 1
    
    ##############################
    
    test_newdata <- data.frame(size = min(training_df$size):max(training_df$size),
                               dummy_var = 0,
                               block = 4,
                               effort1 = training_df$effort1[1],
                               effort2 = training_df$effort2[1])
    
    test_newdata$fit <- predict(gam_logit, newdata = test_newdata, type = "response")
    
    ggplot() +
      geom_path(data = test_newdata,
                mapping = aes(x = size, y = fit))
    
    ##############################
    
    return(validation_df)
  }
  
  doParallel::stopImplicitCluster()
  
  output <- do.call("rbind", cv_results)
  
  output <- output[, -which(names(output) == "dummy_var")]
  
  return(output)
  
  
}
