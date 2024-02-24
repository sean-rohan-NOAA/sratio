#' Simulate catch using a Dirichlet-multinomial distribution
#' 
#' Simulate catch-at-size as the product of selectivity and numbers-at-length using a Dirichlet-multinomial distribution.
#' 
#' @param size Vector of sizes
#' @param abundance Vector of numbers-at-size
#' @param selectivity_opts A vector of selectivity-at-size
#' @param n_size_samples Number of individuals captured per sample
#' @export

catch_at_size <- function(size,
                          abundance,
                          selectivity,
                          n_size_samples) {
  
  abundance_at_size <- data.frame(size = size,
                                  abundance = abundance)
  
  # Draw age/length sample from a multinomial distribution
  # Probability of capture = Dirichlet abundance * selectivity
  comp_df <- stats::rmultinom(1, 
                              size = n_size_samples, 
                              prob = gtools::rdirichlet(n = 1, 
                                                        alpha = abundance_at_size$abundance * selectivity))
  
  sample_comp_df <- abundance_at_size[rep(row.names(abundance_at_size), comp_df), 1:2]
  
  sample_comp_df$specimen_number <- sample(1:nrow(sample_comp_df), 
                                           size = nrow(sample_comp_df), replace = FALSE)
  
  sample_comp_df <- sample_comp_df[, -which(names(sample_comp_df) == "abundance")]
  
  rownames(sample_comp_df) <- 1:nrow(sample_comp_df)
  
  return(sample_comp_df)
  
}