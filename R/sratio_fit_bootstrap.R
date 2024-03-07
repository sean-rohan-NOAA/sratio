#' Bootstrap selectivity ratio
#' 
#' Uses the output of two_stage_bootstrap to fit binomial or beta generalized additive mixed models for catch comparison rate to generate selectivity ratio estimates.
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

sratio_fit_bootstrap <- function(x, treatment_order, size_col, block_col, treatment_col, count_col, effort_col, gam_family, k, scale_method = "sv", n_cores = 1) {
  
  # Make treatment_order a factor and retain the rank order from input
  if(class(treatment_order) != "factor") {
    treatment_order <- factor(treatment_order, levels = treatment_order)
  }
  
  # Function to reformat list of bootstrap samples
  format_sratio_bootstrap <- function(x, size_col, block_col, treatment_col, count_col, effort_col, treatment_order) {
    
    names(x)[match(c(size_col, block_col, treatment_col, count_col, effort_col), table = names(x))] <- c("size", "block", "treatment", "count", "effort")
    
    if(class(x$block) != "factor") {
      x$block <- factor(x$block)
    }
    
    if(class(x$treatment) != "factor") {
      x$treatment <- factor(x$treatment, levels = treatment_order)
    }
    
    lengths <- x |>
      dplyr::mutate(treatment_val = paste0("n_", treatment)) |>
      dplyr::select(block, size, count, treatment_val) |>
      tidyr::pivot_wider(names_from = treatment_val, values_from = count, values_fill = 0) |>
      as.data.frame()
    
    effort <- x |>
      dplyr::mutate(effort_val = paste0("effort_", treatment)) |>
      dplyr::select(effort_val, effort, block) |>
      unique() |>
      tidyr::pivot_wider(names_from = effort_val, values_from = effort, values_fill = 0) |>
      as.data.frame()
    
    combined <- dplyr::inner_join(lengths, effort, by = "block")
    
    combined$p <- sratio::selectivity_ratio(count1 = combined[[paste0("n_", treatment_order[1])]], 
                                            count2 = combined[[paste0("n_", treatment_order[2])]], 
                                            effort1 = combined[[paste0("effort_", treatment_order[1])]], 
                                            effort2 = combined[[paste0("effort_", treatment_order[2])]])$p12
    
    combined$p_scaled <- sratio::scale_for_betareg(combined$p, method = "sv")
    combined$dummy_var <- 1
    
    return(combined)
    
  }
  
  # Format bootstrap data.frames
  x <- lapply(X = x, 
              FUN = format_sratio_bootstrap, 
              size_col = size_col, 
              block_col = block_col, 
              treatment_col = treatment_col, 
              count_col = count_col, 
              effort_col = effort_col,
              treatment_order = treatment_order)
  
  # Get prediction range for sizes
  size_values <- seq(min(unlist(lapply(x, FUN = function(z) {min(z[["size"]])}))), 
                 max(unlist(lapply(x, FUN = function(z) {max(z[["size"]])}))),
                 by = 1)
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  bootstrap_output <- foreach::foreach(iter = 1:length(x),
                                       .packages = c("mgcv", "dplyr")) %dopar% {
  
    boot_df <- x[[iter]]
    
    if(gam_family == "binomial") {
      model <- mgcv::gam(formula = p ~ s(block, bs = 're', by = dummy_var) + s(size, k = k, bs = 'tp'),
                         family = stats::binomial(link = "logit"),
                         data = boot_df)
    }
    
    if(gam_family == "beta") {
    model <- mgcv::gam(formula = p_scaled ~ s(size, k = k, bs = 'tp') + s(block, bs = 're', by = dummy_var),
                       family = mgcv::betar(link = "logit"),
                       data = boot_df)
    }
    
    fit_df <- data.frame(size = size_values,
                         block = x[[iter]]$block[1],
                         dummy_var = 0) # random effects off
    
    # Calculate selectivity ratio
    fit_df$p12 <- predict(model, newdata = fit_df, type = "response")
    fit_df$s21 <- 1/fit_df$p12-1
    
    return(fit_df)
    
  }
  
  doParallel::stopImplicitCluster()
  
  results <- do.call("rbind", bootstrap_output)
  
  names(results)[match(c("size", "block"), table = names(results))] <- c(size_col, block_col)
  
  return(results)
  
}
