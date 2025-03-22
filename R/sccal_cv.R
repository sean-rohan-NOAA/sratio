#' Selectivity conditioned on catch-at-length (SCCAL) cross validation
#' 
#' Block-level cross validation for selectivity condition on catch-at-length Poisson, negative binomial, and Tweedie generalized additive models.
#' 
#' @param gam_family Valid family for `mgcv::gam()`
#' @param gam_formula Valid formula for `mgcv::gam()`
#' @param size1  Numeric vector of sizes or ages for gear #1
#' @param size2  Numeric vector of sizes or ages for gear #2
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param effort1 Numeric vector of effort for for gear #1
#' @param effort2 Numeric vector of effort for for gear #2
#' @param block1 Sample block (i.e. paired sample) for gear #1
#' @param block2 Sample block (i.e. paired sample) for gear #2
#' @param sampling_factor1 Numeric vector of sampling factors for count from gear #1 (i.e., the inverse of the proportion of the catch that was sampled).
#' @param sampling_factor2 Numeric vector of sampling factors for count from gear #2 (i.e., the inverse of the proportion of the catch that was sampled).
#' @param treatment_name1 1L vector for the name of the first treatment
#' @param treatment_name2 1L vector for the name of the second treatment
#' @param k k to use for GAMs. Automatically set to the minimum of 8 or 3 less than the number of unique values in size.
#' @param n_cores Number of cores to use for parallel processing.
#' @returns Returns a data.frame with predicted catch for Poisson, negative binomial, and Tweedie models.
#' @export

sccal_cv <- function(gam_family = poisson(link = "log"),
                     gam_formula = formula(count ~ s(size, k = k, bs = "tp", by = treatment) + s(block, bs = "re") + 
                                             offset(I(log(effort)))),
                     size1, 
                     size2, 
                     count1, 
                     count2, 
                     effort1 = 1, 
                     effort2 = 1, 
                     block1, 
                     block2, 
                     sampling_factor1 = 1, 
                     sampling_factor2 = 1,
                     treatment_name1, 
                     treatment_name2, 
                     k = NULL, 
                     n_cores = 1) {
  
  if(length(sampling_factor1) == 1) {
    sampling_factor1 <- rep(sampling_factor1, length(size1))
  }
  
  if(length(sampling_factor2) == 1) {
    sampling_factor2 <- rep(sampling_factor2, length(size2))
  }
  
  if(length(effort1) == 1) {
    effort1 <- rep(effort1, length(size1))
  }
  
  if(length(effort2) == 1) {
    effort2 <- rep(effort2, length(size2))
  }
  
  model_df <- data.frame(
    count = c(count1 * sampling_factor1, 
              count2 * sampling_factor2),
    sampling_factor = c(sampling_factor1, sampling_factor2),
    size = c(size1, size2),
    effort = c(effort1, effort2),
    block = c(block1, block2),
    treatment = c(rep(treatment_name1, length(count1)), rep(treatment_name2, length(count2)))
  )
  
  if(class(model_df$block) != "factor") {
    model_df$block <- factor(model_df$block)
  }

  if(class(model_df$treatment) != "factor") {
    model_df$treatment <- factor(model_df$treatment)
  }
  
  sccal_fit_gamm <- 
    function(data = data, 
             k = k,
             gam_formula = gam_formula,
             gam_family = gam_family) {
      
      model <- mgcv::gam(formula = gam_formula, family = gam_family, data = data)
      
      return(model)
      
    }
  
  # Setup four clusters and folds for each match-ups
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  index <- as.numeric(interaction(model_df$block, model_df$treatment))
  
  cv_results <- foreach::foreach(fold = 1:max(index), .packages = "mgcv") %dopar% {
    
    training_df <- model_df[-which(index == fold), ]
    validation_df <- model_df[which((index == fold)), ]
    
    mod <- sccal_fit_gamm(data = training_df,
                          k = k,
                          gam_formula = gam_formula,
                          gam_family = gam_family)
    
    validation_df$fit <- predict(object = mod, 
                                 newdata = validation_df, 
                                 type = "response",
                                 exclude = "s(block)")
    
    return(validation_df)
  }
  
  output_cv <- do.call("rbind", cv_results)
  
  doParallel::stopImplicitCluster()
  
  return(
    list(
      cv_results = output_cv,
      model_settings = 
        list(
          gam_family = gam_family,
          k = k
        ),
      data = model_df
    )
  )
  
}
