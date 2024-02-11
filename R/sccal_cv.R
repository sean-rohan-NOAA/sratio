#' Selectivity conditioned on catch-at-length (SCCAL) cross validation
#' 
#' Block-level cross validation for selectivity condition on catch-at-length Poisson, negative binomial, and Tweedie generalized additive models.
#' 
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param effort1 Numeric vector of effort for for gear #1
#' @param effort2 Numeric vector of effort for for gear #2
#' @param size1  Numeric vector of sizes or ages for gear #1
#' @param size2  Numeric vector of sizes or ages for gear #2
#' @param block1 Sample block (i.e. paired sample) for gear #1
#' @param block2 Sample block (i.e. paired sample) for gear #2
#' @param treatment_name1 1L vector for the name of the first treatment
#' @param treatment_name2 1L vector for the name of the second treatment
#' @param k k to use for GAMs. Automatically set to the minimum of 8 or 3 less than the number of unique values in size.
#' @param n_cores Number of cores to use for parallel processing.
#' @returns Returns a data.frame with predicted catch for Poisson, negative binomial, and Tweedie models.
#' @export

sccal_cv <- function(count1, count2, effort1, effort2, size1, size2, block1, block2, treatment_name1, treatment_name2, k = NULL, n_cores = 1) {
  
  model_df <- data.frame(
    count = c(count1, count2),
    size = c(size1, size2),
    effort = c(effort1, effort2),
    block = c(block1, block2),
    treatment = c(rep(treatment_name1, length(count1)), rep(treatment_name2, length(count2))),
    dummy_var = 1
  )
  
  if(class(model_df$block) != "factor") {
    model_df$block <- factor(model_df$block)
  }

  if(class(model_df$treatment) != "factor") {
    model_df$treatment <- factor(model_df$treatment)
  }
  
  # Setup four clusters and folds for each matchups
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  index <- as.numeric(interaction(model_df$block, model_df$treatment))
  
  cv_results <- foreach::foreach(fold = 1:max(index), .packages = "mgcv") %dopar% {
    
    training_df <- model_df[-which(index == fold), ]
    validation_df <- model_df[which((index == fold)), ]
    
    mod_tw <- mgcv::gam(formula = count ~ s(size, k = k, bs = "tp", by = treatment) + 
                          s(block, 
                            by = dummy_var,
                            bs = "re") + 
                          offset(I(log(effort))), 
                        data = training_df,
                        family = tw(link = "log"))
    
    mod_nb <- mgcv::gam(formula = count ~ s(size, k = k, bs = "tp", by = treatment) + 
                          s(block, 
                            by = dummy_var,
                            bs = "re") + 
                          offset(I(log(effort))), 
                        data = training_df,
                        family = nb(link = "log"))
    
    mod_poisson <- mgcv::gam(formula = count ~ s(size, k = k, bs = "tp", by = treatment) + 
                               s(block, 
                                 by = dummy_var,
                                 bs = "re") + 
                               offset(I(log(effort))), 
                             data = training_df,
                             family = poisson(link = "log"))
    
    validation_df$tw <- predict(mod_tw, newdata = validation_df, type = "response")
    validation_df$nb <- predict(mod_nb, newdata = validation_df, type = "response")
    validation_df$poisson <- predict(mod_poisson, newdata = validation_df, type = "response")
    
    return(validation_df)
  }
  
  results <- do.call("rbind", cv_results)
  
  doParallel::stopImplicitCluster()
  
  return(results)
  
}
