#' Two-stage bootstrap draw with resampling
#' 
#' @param x A vector of the variable to be sampled.
#' @param group A vector of grouping variables for x.
#' @param frequency A numeric vector giving the frequency of x. Used for probabalistic sample draws. Default = NULL assigns equal sampling probability to each item x.
#' @param draws Number of group samples to draw
#' @returns A list containing a list of samples and associated group from which the samples were drawn.
#' @export

two_stage_bootstrap <- function(x, group, frequency = NULL, draws) {
  
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


#' Bootstrap draw with resampling
#' 
#' @param x A vector of the variable to be sampled
#' @param group A vector of grouping variables for x.
#' @param frequency A numeric vector giving the frequency of x. Used for probabilistic sample draws. Default = NULL assigns equal sampling probability to each item x.
#' @param draw_group Number of group samples to draw
#' @returns A list containing a list of samples and associated group from which the samples were drawn.
#' @export

nested_bootstrap <- function(x, frequency = NULL, group, draw_group, draw_index) {
  
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
