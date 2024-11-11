#' Selectivity ratio cross-validation
#' 
#' Block-level cross validation for selectivity ratio binomial and beta binomial generalized additive models.
#' 
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param effort1 Numeric vector of effort for for gear #1
#' @param effort2 Numeric vector of effort for for gear #2
#' @param sampling_factor1 Numeric vector of sampling factors for count from gear #1 (i.e., the inverse of the proportion of the catch that was sampled).
#' @param sampling_factor2 Numeric vector of sampling factors for count from gear #2 (i.e., the inverse of the proportion of the catch that was sampled).
#' @param obs_weight_options A list indicating which method to use to weight observations. The default is to not weight observations (obs_weight_options = NULL). The current alternative is to weight observations by the inverse of the residual variance based on total count, resid ~ s(count1 + count 2) (obs_weight_options = "count").
#' @param size  Numeric vector of sizes or ages
#' @param block Sample block (i.e. paired sample)
#' @param k k to use for GAMs. Automatically set to the minimum of 8 or 3 less than the number of unique values in size.
#' @param scale_method Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg
#' @param n_cores Number of cores to use for parallel processing.
#' @export

sratio_cv <- function(size, 
                      count1, 
                      count2, 
                      effort1, 
                      effort2, 
                      sampling_factor1 = 1, 
                      sampling_factor2 = 1, 
                      obs_weight_options = NULL,
                      block, 
                      k = NULL, 
                      scale_method = "sv", 
                      n_cores = 1) {

  unique_blocks <- unique(block)
  
  p <- numeric(length = length(count1))
  s <- numeric(length = length(count1))
  
  if(length(sampling_factor1) == 1) {
    sampling_factor1 <- rep(sampling_factor1, length(count1))
  }
  
  if(length(sampling_factor2) == 1) {
    sampling_factor2 <- rep(sampling_factor2, length(count2))
  }
  
  for(jj in 1:length(unique_blocks)) {
    
    s_ratio_df <- suppressMessages(
      selectivity_ratio(size = size[block == unique_blocks[jj]],
                        count1 = count1[block == unique_blocks[jj]], 
                        count2 = count2[block == unique_blocks[jj]],
                        effort1 = effort1[block == unique_blocks[jj]], 
                        effort2 = effort2[block == unique_blocks[jj]],
                        sampling_factor1 = sampling_factor1[block == unique_blocks[jj]],
                        sampling_factor2 = sampling_factor2[block == unique_blocks[jj]]
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
  sampling_factor1 <- sampling_factor1[!is.na(p)]
  sampling_factor2 <- sampling_factor2[!is.na(p)]
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
                         total_count = count1 + count2,
                         effort1 = effort1,
                         effort2 = effort2,
                         sampling_factor1 = sampling_factor1,
                         sampling_factor2 = sampling_factor2,
                         p = p,
                         p_scaled = p_scaled,
                         s = s,
                         dummy_var = 1)
  
  # Setup four clusters and folds for each block
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  folds <- caret::groupKFold(group = block)
  
  cv_results <- foreach::foreach(fold = folds, .packages = c("mgcv", "dplyr")) %dopar% {
    
    training_df <- model_df[fold, ]
    validation_df <- model_df[-fold, ]
    # validation_df$dummy_var <- 0
    
    # Add in dummy station variable for predictions, to be added back in for output
    out_matchup <- validation_df$block[1]
    validation_df$block <- training_df$block[1]
    
    # No weights
    weight_logit <- weight_beta <- NULL
    
    # Weight observations by count, based on either the relationship between total count and residuals or just by count
    # if(is(obs_weight_options, "list")) {
    
    if(is.numeric(obs_weight_options$max_count)) {
      
      training_df$total_count <- 
        dplyr::if_else(
          training_df$total_count > obs_weight_options$max_count,
          obs_weight_options$max_count,
          training_df$total_count
        )
      
    }
    
    if(!is.null(obs_weight_options$method)) {
      
      if(obs_weight_options$method == "count") { # By count
        
        weight_logit <- weight_beta <- training_df$total_count
        
      }
      
      if(obs_weight_options$method == "residuals_by_count") { # By residuals
        
        # Calculate absolute or squared residuals
        resid_fun <- switch(
          obs_weight_options$residual_type,
          "absolute" = function(x) abs(resid(x)),
          "squared" = function(x) resid(x)^2)
        
        training_df$transformed_resid_logit <- resid_fun(gam_logit)
        
        training_df$transformed_resid_beta <- resid_fun(gam_beta)
        
        # Fit model model residual and predict
        gam_logit_dev_mod <- mgcv::gam(transformed_resid_logit ~ s(total_count, bs = "tp"), data = training_df)
        
        gam_beta_dev_mod <- mgcv::gam(transformed_resid_beta ~ s(total_count, bs = "tp"), data = training_df)
        
        # Estimate residual fit for each observation
        training_df$fit_inv_var_logit <- 1/predict(object = gam_logit_dev_mod, 
                                                   newdata = training_df,
                                                   type = "response")
        
        training_df$fit_inv_var_beta <- 1/predict(object = gam_beta_dev_mod,
                                                  newdata = training_df,
                                                  type = "response")
        
        # Transform estimate of SD to variance when model is fit to absolute residuals  
        if(obs_weight_options$residual_type == "absolute") {
          
          training_df$fit_inv_var_logit <- training_df$fit_inv_var_logit^2
          
          training_df$fit_inv_var_beta <- training_df$fit_inv_var_beta^2
          
        }
        
        # Calculate normalized weights for each observation
        weight_logit <- training_df$fit_inv_var_logit #/ mean(training_df$fit_inv_var_logit)
        
        weight_beta <- training_df$fit_inv_var_beta #/ mean(training_df$fit_inv_var_beta)
        
      }
    }
      
    # }
    
    gam_logit <- mgcv::gam(p ~ s(size, bs = "tp", k = k) + s(block, bs = "re"),
                           data = training_df,
                           weights = weight_logit,
                           family = binomial(link = "logit"))
    
    gam_beta <- mgcv::gam(p_scaled ~ s(size, bs = "tp", k = k) + s(block, bs = "re"),
                          weights = weight_beta,
                          data = training_df,
                          family = mgcv::betar(link = "logit"))
    
    fitted_logit <- predict(gam_logit, newdata = validation_df, type = "response", exclude = "s(size)")
    
    fitted_beta <- predict(gam_beta, newdata = validation_df, type = "response", exclude = "s(size)")
    
    validation_df$cv_fit_logit <- fitted_logit
    validation_df$cv_fit_beta <- fitted_beta
    
    # Reset match-ups for fitting final models
    validation_df$block <- out_matchup
    
    return(validation_df)
  }
  
  doParallel::stopImplicitCluster()
  
  output <- do.call("rbind", cv_results)
  
  return(output)
  
}
