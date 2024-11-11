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
#' @param sampling_factor_col Character vector; name of the column containing sampling factor (count multiplier, where estimated catch-at-length = count * sampling_factor)
#' @param obs_weight_control A list indicating which method to use to weight observations. See help documentation `sratio_fit_gamm()` for information about options (?sratio_fit_gamm) .
#' @param gam_family "binomial" or "beta"
#' @param k k for mgcv spline for size
#' @param scale_method Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg
#' @param n_cores Number of cores to use for parallel processing.
#' @export

sratio_fit_bootstrap <- function(x, 
                                 treatment_order, 
                                 size_col, 
                                 block_col, 
                                 treatment_col, 
                                 count_col, 
                                 effort_col, 
                                 sampling_factor_col = NULL,
                                 obs_weight_control = list(method = "none", 
                                                           max_count = Inf,
                                                           residual_type = NA,
                                                           normalize_weights = FALSE),
                                 gam_family, 
                                 k, 
                                 scale_method = "sv", 
                                 n_cores = 1) {
  
  # Make treatment_order a factor and retain the rank order from input
  if(class(treatment_order) != "factor") {
    treatment_order <- factor(treatment_order, levels = treatment_order)
  }
  
  # Function to reformat data.frames within a list containing bootstrap samples
  format_sratio_bootstrap <- function(x, 
                                      size_col, 
                                      block_col, 
                                      treatment_col, 
                                      count_col, 
                                      effort_col, 
                                      sampling_factor_col, 
                                      treatment_order,
                                      scale_method) {
    
    names(x)[
      match(
        c(size_col, block_col, treatment_col, count_col, effort_col, sampling_factor_col), 
        table = names(x)
      )
    ] <- c("size", "block", "treatment", "count", "effort", "sampling_factor")
    
    if(!is(x$block, "factor")) {
      x$block <- factor(x$block)
    }
    
    if(!is(x$treatment, "factor")) {
      x$treatment <- factor(x$treatment, levels = treatment_order)
    }
    
    # Append prefixes to treatment levels
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
    
    sampling_factor <- x |>
      dplyr::mutate(sampling_factor_val = paste0("sampling_factor_", treatment)) |>
      dplyr::select(sampling_factor_val, size, sampling_factor, block) |>
      tidyr::pivot_wider(names_from = sampling_factor_val, values_from = sampling_factor, values_fill = 1) |>
      unique() |>
      as.data.frame()
    
    combined <- dplyr::inner_join(lengths, effort, by = "block") |>
      dplyr:::inner_join(sampling_factor, by = c("block", "size"))
    
    sr_values <- sratio::selectivity_ratio(count1 = combined[[paste0("n_", treatment_order[1])]], 
                                           count2 = combined[[paste0("n_", treatment_order[2])]], 
                                           effort1 = combined[[paste0("effort_", treatment_order[1])]], 
                                           effort2 = combined[[paste0("effort_", treatment_order[2])]],
                                           sampling_factor1 = combined[[paste0("sampling_factor_", treatment_order[1])]],
                                           sampling_factor2 = combined[[paste0("sampling_factor_", treatment_order[2])]])
    
    combined$total_count <- sr_values$count1 + sr_values$count2
    
    combined$p <- sr_values$p12
    
    combined$p_scaled <- sratio::scale_for_betareg(combined$p, method = scale_method)
    
    return(combined)
    
  }
  
  if(is.null(sampling_factor_col)) {
    sampling_factor_col <- "sampling_factor"
    x$sampling_factor <- 1
  }
  
  # Format bootstrap data.frames
  x <- lapply(X = x, 
              FUN = format_sratio_bootstrap, 
              size_col = size_col, 
              block_col = block_col, 
              treatment_col = treatment_col, 
              count_col = count_col, 
              effort_col = effort_col,
              sampling_factor_col = sampling_factor_col,
              treatment_order = treatment_order,
              scale_method = scale_method)
  
  # Get prediction range for sizes
  size_values <- seq(min(unlist(lapply(x, 
                                       FUN = function(z) {floor(min(z[["size"]]))}))), 
                 max(unlist(lapply(x, 
                                   FUN = function(z) {ceiling(max(z[["size"]]))}))),
                 by = 1)
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  bootstrap_output <- foreach::foreach(iter = 1:length(x),
                                       .packages = c("mgcv", "dplyr")) %dopar% {
  
    boot_df <- x[[iter]]
    
    # Set model family based on user-specified settings
    if(gam_family == "binomial") {
      
      fit_family <- binomial(link = "logit")
      
    }
    
    if(gam_family == "beta") {
      
      fit_family <- binomial(link = "logit")
      
    }
    
    # Fit model
    model <- 
      sratio_fit_gamm(data = boot_df,
                             k = k,
                             gam_formula = 
                               p ~ s(size, bs = "tp", k = k) + s(block, bs = "re"),
                             gam_family = fit_family, 
                             obs_weight_control = obs_weight_control)$mod
    
    # Create prediction data.frame. Includes unused block value to avoid an error; doesn't get used with random effects off
    fit_df <- data.frame(size = size_values,
                         block = -999) 
    
    # Calculate selectivity ratio
    fit_df$p12 <- predict(model, 
                          newdata = fit_df,
                          exclude = "s(block)", # random effect of block is off
                          type = "response")
    fit_df$s12 <- fit_df$p12/(1-fit_df$p12)
    fit_df$draw <- iter
    
    return(fit_df)
    
  }
  
  doParallel::stopImplicitCluster()
  
  results <- do.call("rbind", bootstrap_output)
  
  names(results)[match(c("size", "block"), table = names(results))] <- c(size_col, block_col)
  
  return(results)
  
}
