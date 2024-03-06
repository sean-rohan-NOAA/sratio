#' Two-stage bootstrap from paired samples
#' 
#' Two-stage bootstrap resampling (by block/pair then individual within each pair), where samples are drawn from paired treatments.
#' 
#' @param count1 Numeric vector of catch-at-size/age for gear #1
#' @param count2 Numeric vector of catch-at-size/age for gear #2
#' @param size1  Numeric vector of sizes or ages for gear #1
#' @param size2  Numeric vector of sizes or ages for gear #2
#' @param block1 Sample block for gear #1
#' @param block2 Sample block for gear #2
#' @param n_draws Number of samples to draw.
#' @param seed random number generator seed
#' @param treatment_name1 Name of the treatment for the first gear as a number or character.
#' @param treatment_name2 Name of the treatment for the second gear as a number or character.
#' @export

two_stage_bootstrap <- function(count1,
                                count2,
                                size1,
                                size2,
                                block1,
                                block2,
                                n_draws,
                                seed = NULL,
                                treatment_name1 = 1,
                                treatment_name2 = 2) {
  
  set.seed(seed)
  boot_1 <- boot_stage_one(x = size1,
                           group = block1, 
                           frequency = count1,
                           draws = n_draws)
  
  set.seed(seed)
  boot_2 <- boot_stage_two(x = size2, 
                           frequency = count2, 
                           group = block2, 
                           draw_group = boot_1$draws$group, 
                           draw_index = boot_1$draws$draw)
  
  boot_samples <- vector(mode = "list", length = n_draws)
  
  for(kk in 1:n_draws) {
    
    boot_samples[[kk]] <- rbind(
      dplyr::mutate(boot_1$samples[[kk]],
                    treatment = treatment_name1) |>
        dplyr::rename(size = x, 
                      block = group),
      dplyr::mutate(boot_2$samples[[kk]],
                    treatment = treatment_name2) |>
        dplyr::rename(size = x, 
                      block = group)
    ) |>
      dplyr::group_by(block, treatment, size) |>
      dplyr::summarize(count = n(), .groups = "keep") |>
    as.data.frame()
    
  }
  
  return(boot_samples)
  
}


#' First stage of two-stage bootstrap
#' 
#' Draw from paired treatment blocks and samples from first treatment level
#' 
#' @param x A vector of the variable to be sampled.
#' @param group A vector of grouping variables for x.
#' @param frequency A numeric vector giving the frequency of x. Used for probabalistic sample draws. Default = NULL assigns equal sampling probability to each item x.
#' @param draws Number of group samples to draw
#' @returns A list containing a list of samples and associated group from which the samples were drawn.
#' @export

boot_stage_one <- function(x, group, frequency = NULL, draws) {
  
  if(is.null(frequency)) {
    frequency <- rep(1, length(x))
  }
  
  n_group <- length(unique(group))
  
  samples_list <- vector(mode = "list", length = draws)
  
  group_draw <- vector()
  iter <- numeric()
  
  for(jj in 1:draws) {
    
    boot_group <- sort(sample(x = unique(group),
                              size = n_group,
                              replace = TRUE))
    
    group_draw <- c(group_draw, boot_group)
    iter <- c(iter, rep(jj, n_group))
    
    group_var <- "Go Huskies"
    
    samples_df <- data.frame()
    
    for(ii in 1:n_group) {
      
      if(boot_group[ii] != group_var) {
        group_var <- boot_group[ii]
        sel_x <- x[group == group_var]
        sel_freq <- frequency[group == group_var]
        
      }
      
      sel_n <- round(sum(sel_freq))
      
      if(length(sel_x) == 1) {
        draw_df <- data.frame(x = sel_x,
                              group = group_var)
      } else {
        draw_df <- data.frame(x = sample(sel_x, 
                                         size = sel_n, 
                                         replace = TRUE, 
                                         prob = sel_freq),
                              group = group_var)
      }
      
      samples_df <- rbind(samples_df, draw_df)
      
    }
    
    samples_list[[jj]] <- samples_df
    
  }
  
  return(list(samples = samples_list, draws = data.frame(draw = iter, group = group_draw)))
  
}


#' Second stage of two-stage bootstrap
#' 
#' Draw from second treatment level using blocks selected for the first stage of two-stage bootstrap
#' 
#' @param x A vector of the variable to be sampled
#' @param group A vector of grouping variables for x.
#' @param frequency A numeric vector giving the frequency of x. Used for probabilistic sample draws. Default = NULL assigns equal sampling probability to each item x.
#' @param draw_group Number of group samples to draw
#' @param draw_index Index value for bootstrap draw (i.e., samples drawn the first sample = 1, second sample = 2, etc.)
#' @returns A list containing a list of samples and associated group from which the samples were drawn.
#' @export

boot_stage_two <- function(x, frequency = NULL, group, draw_group, draw_index) {
  
  if(is.null(frequency)) {
    frequency <- rep(1, length(x))
  }
  
  n_draws <- length(unique(draw_index))
  
  samples_list <- vector(mode = "list", length = n_draws)
  
  group_var <- "Go Huskies"
  
  for(ii in 1:length(draw_group)) {
    
    if(draw_group[ii] != group_var) {
      group_var <- draw_group[ii]
      sel_x <- x[group == group_var]
      sel_freq <- frequency[group == group_var]
      
    }
    
    sel_n <- round(sum(sel_freq))
    
    if(length(sel_x) == 1) {
      
      samples_list[[draw_index[ii]]] <- rbind(samples_list[[draw_index[ii]]], 
                                              data.frame(x = sel_x,
                                                         group = group_var))
    } else{
      samples_list[[draw_index[ii]]] <- rbind(samples_list[[draw_index[ii]]], 
                                              data.frame(x = sample(sel_x,
                                                                    size = sel_n, 
                                                                    replace = TRUE, 
                                                                    prob = sel_freq),
                                                         group = group_var)
      )
    }
    
  }
  
  return(list(samples = samples_list, 
              group = data.frame(draw = draw_index,
                                 group = draw_group)))
  
}


