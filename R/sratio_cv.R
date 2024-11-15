#' Selectivity ratio cross-validation
#' 
#' Block-level cross validation for selectivity ratio binomial and beta binomial generalized additive models.
#' 
#' @param model_type Model type. Options are "binomial" or "beta"
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param effort1 Numeric vector of effort for for gear #1
#' @param effort2 Numeric vector of effort for for gear #2
#' @param sampling_factor1 Numeric vector of sampling factors for count from gear #1 (i.e., the inverse of the proportion of the catch that was sampled).
#' @param sampling_factor2 Numeric vector of sampling factors for count from gear #2 (i.e., the inverse of the proportion of the catch that was sampled).
#' @param obs_weight_control A list indicating which method to use to weight observations. See help documentation `sratio_fit_gamm()` for information about options (?sratio_fit_gamm) .
#' @param size  Numeric vector of sizes or ages
#' @param block Sample block (i.e. paired sample)
#' @param k k to use for GAMs. Automatically set to the minimum of 8 or 3 less than the number of unique values in size.
#' @param scale_method Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg
#' @param sratio_type Which selectivity ratio calculation should be used? Absolute ("absolute") or relative ("relative")?
#' @param n_cores Number of cores to use for parallel processing.
#' @export

sratio_cv <- function(model_type = "binomial",
                      size, 
                      count1, 
                      count2, 
                      effort1, 
                      effort2, 
                      sampling_factor1 = 1, 
                      sampling_factor2 = 1, 
                      obs_weight_control = 
                        list(method = "none",
                             max_count = Inf,
                             residual_type = NA,
                             normalize_weight = FALSE),
                      block, 
                      k = NULL, 
                      scale_method = "sv", 
                      sratio_type = "absolute",
                      n_cores = 1) {
  
  # Initialize vectors for catch comparison rate (p) and selectivty ratio (s)
  p <- numeric(length = length(count1))
  s <- numeric(length = length(count1))
  
  # Create sampling factor vectors if provided as a 1L numeric
  if(length(sampling_factor1) == 1) {
    sampling_factor1 <- rep(sampling_factor1, length(count1))
  }
  
  if(length(sampling_factor2) == 1) {
    sampling_factor2 <- rep(sampling_factor2, length(count2))
  }
  
  # Set k if not provided
  if(is.null(k)) {
    k <- min(c(8, (length(unique(size))-3)))
  }
  
  # Set model family
  if(model_type == "binomial") {
    gam_family <- binomial(link = "logit")
  }
  
  if(model_type == "beta") {
    gam_family <- betar(link = "logit")
    # Scale p to fall within (the interval supported by (0,1] for beta regression
    p <- scale_for_betareg(p, method = scale_method)
  }
  
  
  # Estimate selectivity ratios for each size within each block
  unique_blocks <- unique(block)
  
  for(jj in 1:length(unique_blocks)) {
    
    s_ratio_df <- suppressMessages(
      selectivity_ratio(
        size = size[block == unique_blocks[jj]],
        count1 = count1[block == unique_blocks[jj]], 
        count2 = count2[block == unique_blocks[jj]],
        effort1 = effort1[block == unique_blocks[jj]], 
        effort2 = effort2[block == unique_blocks[jj]],
        sampling_factor1 = sampling_factor1[block == unique_blocks[jj]],
        sampling_factor2 = sampling_factor2[block == unique_blocks[jj]],
        sratio_type = sratio_type
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
  s <- s[!is.na(p)]
  p <- p[!is.na(p)]

  # Setup data.frame for leave-one-out cross-validation (LOOCV)
  model_df <- data.frame(
    block = factor(block),
    size = size,
    count1 = count1,
    count2 = count2,
    total_count = count1 + count2,
    effort1 = effort1,
    effort2 = effort2,
    sampling_factor1 = sampling_factor1,
    sampling_factor2 = sampling_factor2,
    p = p,
    s = s
  )
  
  # Setup parallelization and folds for LOOCV
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  folds <- caret::groupKFold(group = block)
  
  # Conduct leave one out cross validation
  cv_results <- foreach::foreach(fold = folds, .packages = c("mgcv", "dplyr")) %dopar% {
    
    training_df <- model_df[fold, ]
    validation_df <- model_df[-fold, ]
    
    mod <- sratio_fit_gamm(
      data = training_df,
      k = k,
      gam_formula = p ~ s(size, bs = "tp", k = k) + s(block, bs = "re"),
      gam_family = gam_family, 
      obs_weight_control = obs_weight_control
    )
    
    # Fit without random effect of block
    validation_df$p_fit <- predict(mod$mod, 
                                   newdata = validation_df, 
                                   type = "response", 
                                   exclude = "s(block)")
    validation_df$s_fit <- validation_df$p_fit/(1-validation_df$p_fit)
    
    # Reset match-ups for fitting final models
    
    return(validation_df)
  }
  
  doParallel::stopImplicitCluster()
  
  # Combine individual outputs
  output_cv <- do.call("rbind", cv_results)

  return(
    list(
      cv_results = output_cv,
      model_settings = 
        list(
          model_type = model_type,
          k = k,
          obs_weight_control = obs_weight_control
        ),
      data = model_df
    )
  )
  
}
