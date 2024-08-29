#' Selectivity ratio cross-validation
#' 
#' Block-level cross validation for selectivity ratio binomial and beta binomial models where size is a factor.
#' 
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param effort1 Numeric vector of effort for for gear #1
#' @param effort2 Numeric vector of effort for for gear #2
#' @param size  Numeric vector of sizes or ages
#' @param block Sample block (i.e. paired sample)
#' @param scale_method Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg
#' @param n_cores Number of cores to use for parallel processing.
#' @export

discrete_cv <- function(count1, count2, effort1, effort2, size, block, scale_method = "sv", n_cores = 1) {
  
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
  
  model_df <- data.frame(block = factor(block),
                         size =  factor(size),
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
  
  cv_results <- foreach::foreach(fold = folds, .packages = c("dplyr")) %dopar% {
    
    training_df <- model_df[fold, ]
    validation_df <- model_df[-fold, ]
    validation_df$dummy_var <- 0
    
    validation_df <- validation_df[validation_df$size %in% training_df$size, ]
    
    # Add in dummy station variable for predictions, to be added back in for output
    out_matchup <- validation_df$block[1]
    validation_df$block <- training_df$block[1]
    
    gam_logit <- mgcv::gam(p ~ size + 0 + s(block, bs = "re", by = dummy_var),
                           data = training_df,
                           family = binomial(link = "logit"))
    
    gam_beta <- mgcv::gam(p_scaled ~ size + 0 + s(block, bs = "re", by = dummy_var),
                          data = training_df,
                          family = mgcv::betar(link = "logit"))
    
    fitted_logit <- predict(gam_logit, newdata = validation_df, type = "response")
    fitted_beta <- predict(gam_beta, newdata = validation_df, type = "response")
    
    validation_df$cv_fit_logit <- fitted_logit
    validation_df$cv_fit_beta <- fitted_beta
    
    # Reset matchup and dummy variable for fitting final models
    validation_df$block <- out_matchup
    validation_df$dummy_var <- 1
    
    return(validation_df)
  }
  
  doParallel::stopImplicitCluster()
  
  output <- do.call("rbind", cv_results)
  
  output <- output[, -which(names(output) == "dummy_var")]
  
  return(output)
  
}
