#' Bootstrap selectivity condition on catch-at-length
#' 
#' Uses the output of two_stage_bootstrap to Poisson, negative binomial, or Tweedie generalized additive mixed models to catch at length data, then calculates selectivity ratio.
#' 
#' @param x list containing bootstrap samples, output from two_stage_bootstrap
#' @param treatment_order rank order of treatment levels, must match treatment names in treatment col
#' @param size_col Character vector; name of the size column
#' @param block_col Character vector; name of the sampling block column
#' @param treatment_col Character vector; name of the column containing treatment names
#' @param count_col Character vector; name of the column containing counts
#' @param effort_col Character vector; name of the column containing effort
#' @param gam_family "binomial" or "beta"
#' @param k k for mgcv spline for size
#' @param scale_method Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg
#' @param n_cores Number of cores to use for parallel processing.
#' @export

sccal_fit_bootstrap <- function(x, treatment_order, size_col, block_col, treatment_col, count_col, effort_col, gam_family, k, n_cores = 1) {
  
  default_effort <- 1 # Default area swept value for generating predictions, 1 km2 (arbitrary)
  
  # Make treatment_order a factor and retain the rank order from input
  if(class(treatment_order) != "factor") {
    factor(treatment_order, levels = treatment_order)
  }
  
  format_sccal_bootstrap <- function(x, size_col, block_col, treatment_col, count_col, effort_col, treatment_order) {
    
    names(x)[match(c(size_col, block_col, treatment_col, count_col, effort_col), table = names(x))] <- c("size", "block", "treatment", "count", "effort")
    
    if(class(x$block) != "factor") {
      x$block <- factor(x$block)
    }
    
    if(class(x$treatment) != "factor") {
      x$treatment <- factor(x$treatment, levels = treatment_order)
    }
    
    x$dummy_var <- 1
    
    return(x)
    
  }
  
  # Format bootstrap data.frame for model fitting
  x <- lapply(X = x, 
              FUN = format_sccal_bootstrap,
              size_col = size_col,
              block_col = block_col,
              treatment_col = treatment_col,
              count_col = count_col,
              effort_col = effort_col,
              treatment_order = treatment_order)
  
  
  size_values <- seq(min(unlist(lapply(x, FUN = function(z) {min(z$size)}))), 
                     max(unlist(lapply(x, FUN = function(z) {max(z$size)}))),
                     by = 1)
  
  sel_family <- switch(gam_family,
                       "tw" = tw(link = "log"),
                       "nb" = nb(link = "log"),
                       "poisson" = poisson(link = "log"))
  
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  bootstrap_output <- foreach::foreach(iter = 1:length(x), 
                                       .packages = c("mgcv", "tidyr")) %dopar% {
    
    model <- mgcv::gam(formula = count ~ s(size, bs = "tp", k = k, by = treatment) + s(block, bs = "re", by = dummy_var) + offset(I(log(effort))),
                       family = sel_family,
                       data = x[[iter]])
    
    fit_df <- expand.grid(size = size_values,
                          block = x[[iter]]$block[1],
                          treatment = treatment_order,
                          effort = default_effort, 
                          dummy_var = 0, # random effects off
                          draw = iter)
    
    fit_df$fit <- predict(model, newdata = fit_df, type = "response")
    
    # Calculate selectivity ratio
    fit_df <- fit_df |>
      tidyr::pivot_wider(id_cols = c("size", "block", "draw", "effort"), 
                         names_from = "treatment", 
                         values_from = "fit")
    
    selectivity_results <- sratio::selectivity_ratio(count1 = fit_df[[which(names(fit_df) == treatment_order[1])]],
                                                     count2 = fit_df[[which(names(fit_df) == treatment_order[2])]],
                                                     effort1 = default_effort,
                                                     effort2 = default_effort)
    
    fit_df$p12 <- selectivity_results$p12
    fit_df$s21 <- selectivity_results$s21
    
    return(fit_df)
    
    
  }
  
  doParallel::stopImplicitCluster()
  
  results <- do.call("rbind", bootstrap_output)
  
  names(results)[match(c("size", "block", "effort"), table = names(results))] <- c(size_col, block_col, effort_col)
  
  return(results)
  
}
