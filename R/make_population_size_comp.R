#' Simulate population size composition from mean size at age
#' 
#' Simulate size composition of the population based on numbers-at-age, size-at-age, and standard deviation in size-at-age.
#' 
#' @param age Vector of age
#' @param size Vector of mean length at age
#' @param abundance Vector of abundances at size
#' @param sd_at_age_df Data frame containing standard deviation at age 1 and 20 for the chosen ensemble.
#' @param size_range Numeric vector of 2L. Maximum and minimum lengths
#' @export

make_population_size_comp <- function(age,
                                 size,
                                 abundance,
                                 sd_at_age_df,
                                 size_range = c(3, 127)) {
  
  # Linear approximation of standard deviation at age
  sd_at_age_vec <- sd_at_age_df$sd_at_age1 + (((age-1) * (sd_at_age_df$sd_at_age20-sd_at_age_df$sd_at_age1))/(20-1))
  
  # Standard deviation at age zero
  sd_at_age_vec[1] <- 0.5
  
  for(i in 1:length(abundance)) {
    # Handling abundance < 1,000,000 fish
    if(abundance[i] < 1000) {
      size <- round((rnorm(abundance[i]*1000, size[i],sd_at_age_vec[i])))
    } else {
      size <- round((rnorm(abundance[i], size[i],sd_at_age_vec[i])))
    }
    
    # Maximum length at 127 cm
    while(any(size > size_range[2])) {
      print(paste0(size[i], " ", length(size[size > size_range[2]])))
      size[size > size_range[2]] <- round(rnorm(length(size[size > size_range[2]]), size[i],sd_at_age_vec[i]))
    }
    
    # Minimum length at 3 cm
    while(any(size < size_range[1])) {
      print(paste0(size[i], " ", length(size[size < size_range[1]])))
      size[size < size_range[1]] <- round(rnorm(length(size[size < size_range[1]]), size[i],sd_at_age_vec[i]))
    }
    
    out_df <- as.data.frame(table(size)) %>%
      dplyr::mutate(age = age[i])
    if(abundance[i] < 1000) {
      out_df$Freq <- out_df$Freq/1000
    }
    
    if(i == 1) {
      size_df <- out_df
    } else {
      size_df <- dplyr::bind_rows(size_df, out_df)
    }
  }
  
  size_wide_df <- size_df %>% 
    tidyr::pivot_wider(id_cols = "size", 
                       names_from = "age", 
                       values_from = "Freq", 
                       values_fill = 0)
  
  size_sums <- data.frame(size_adj = size_wide_df$size,
                        abundance = rowSums(size_wide_df[,2:ncol(size_wide_df)]))
  
  return(size_sums)
  
}