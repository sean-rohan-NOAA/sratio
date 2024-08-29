#' Bootstrap selectivity ratio
#' 
#' Uses the output of two_stage_bootstrap to fit discrete size/age binomial or beta generalized additive mixed models for catch comparison rate to generate selectivity ratio estimates.
#' 
#' @param x list containing bootstrap samples, output from two_stage_bootstrap
#' @param treatment_order rank order of treatment levels, must match treatment names in treatment col
#' @param size_col Character vector; name of the size/age column
#' @param block_col Character vector; name of the sampling block column
#' @param treatment_col Character vector; name of the column containing treatment names
#' @param count_col Character vector; name of the column containing counts
#' @param effort_col Character vector; name of the column containing effort
#' @param gam_family "binomial" or "beta"
#' @param scale_method Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg
#' @param n_cores Number of cores to use for parallel processing.
#' @export

discrete_fit_bootstrap <- function(x, treatment_order, size_col, block_col, treatment_col, count_col, effort_col, gam_family, scale_method = "sv", n_cores = 1) {
  
  # Make treatment_order a factor and retain the rank order from input
  if(class(treatment_order) != "factor") {
    treatment_order <- factor(treatment_order, levels = treatment_order)
  }
  
  # Function to reformat list of bootstrap samples
  format_discrete_bootstrap <- function(x, size_col, block_col, treatment_col, count_col, effort_col, treatment_order) {
    
    names(x)[match(c(size_col, block_col, treatment_col, count_col, effort_col), table = names(x))] <- c("size", "block", "treatment", "count", "effort")
    
    if(class(x$block) != "factor") {
      x$block <- factor(x$block)
    }
    
    if(class(x$treatment) != "factor") {
      x$treatment <- factor(x$treatment, levels = treatment_order)
    }
    
    if(class(x$size) != "factor") {
      x$size <- factor(x$size)
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
      tidyr::pivot_wider(names_from = effort_val, 
                         values_from = effort, 
                         values_fill = 0) |>
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
              FUN = format_discrete_bootstrap, 
              size_col = size_col, 
              block_col = block_col, 
              treatment_col = treatment_col, 
              count_col = count_col, 
              effort_col = effort_col,
              treatment_order = treatment_order)
  
  # Get prediction range for sizes
  size_values <- seq(
    min(
      unlist(
        lapply(x, 
               FUN = function(z) {
                 min(
                   as.numeric(
                     as.character(
                       z[["size"]]
                     )
                   )
                 )
               }
        )
      )
    ), 
    max(
      unlist(
        lapply(x, 
               FUN = function(z) {
                 max(
                   as.numeric(
                     as.character(
                       z[["size"]]
                     )
                   )
                 )
               }
        )
      )
    ),
    by = 1)
  
  size_bins <- sort(unique(unlist(lapply(x, FUN = function(z) {unique(z[["size"]])}))))
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  bootstrap_output <- foreach::foreach(
    iter = 1:length(x),
    .packages = c("mgcv", "dplyr")) %dopar% {
      
      boot_df <- x[[iter]]
      
      if(gam_family == "binomial") {
        model <- mgcv::gam(formula = p ~ size + 0 + s(block, bs = 're', by = dummy_var),
                           family = stats::quasibinomial(link = "logit"),
                           data = boot_df)
      }
      
      if(gam_family == "beta") {
        model <- mgcv::gam(formula = p_scaled ~ size + 0 + s(block, bs = 're', by = dummy_var),
                           family = mgcv::betar(link = "logit"),
                           data = boot_df)
      }
      
      factor_df <- data.frame(size = size_bins,
                              block = x[[iter]]$block[1],
                              dummy_var = 0) |> # random effects off
        dplyr::filter(size %in% unique(boot_df$size))
      
      # Calculate selectivity ratio
      factor_df$p12 <- predict(model, 
                               newdata = factor_df, 
                               type = "response")
      
      factor_df$s12 <- factor_df$p12/(1-factor_df$p12)
      
      fit_df <- data.frame(size = size_values,
                           block = x[[iter]]$block[1],
                           dummy_var = 0)
      
      fit_df$p12 <- approx(x = as.numeric(as.character(factor_df$size)), 
                           y = factor_df$p12,
                           xout = fit_df$size)$y
      
      fit_df$s12 <- fit_df$p12/(1-fit_df$p12)
  
  return(fit_df)
  
}
  
  doParallel::stopImplicitCluster()
  
  results <- do.call("rbind", bootstrap_output)
  
  names(results)[match(c("size", "block"), table = names(results))] <- c(size_col, block_col)
  
  return(results)
  
}
