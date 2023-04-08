#' Run bootstrap sampling for N iterations
#' 
#' @param dat  data.frame containing a columns named p_fishery_adj, MATCHUP, LEN_MIDPOINT, dummy.
#' @param lengths numeric vector of fork lengths for prediction
#' @param iterations Number of iterations to use.
#' @param clusters Number of clusters for parallel computing.
#' @param seed Seed number to use for parallel cluster
#' @param k gam max knots
#' @export

run_bootstrap <- function(dat, lengths, iterations = 1000, clusters = 2, seed, k = 10) {
  
  draw_bootstrap_samples <- function(x) {
    unique_matchup <- unique(x$MATCHUP)
    n_draws <- length(unique_matchup)
    boot_matchups <- sample(x = unique_matchup,
                            size = n_draws,
                            replace = TRUE)
    
    draw_df <- data.frame()
    
    for(ii in 1:n_draws) {
      
      boot_select <- dplyr::filter(x, MATCHUP == boot_matchups[ii])
      
      boot_lengths <- sample(x = 1:nrow(boot_select),
                             size = nrow(boot_select),
                             replace = TRUE)
      
      draw_df <- dplyr::bind_rows(draw_df,
                                  boot_select[boot_lengths, ]
                                  
      )
    }
    return(draw_df)
  }
  
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  parallel::clusterSetRNGStream(cl = cl, iseed = seed)
  
  bootstrap_output <- foreach::foreach(n_iter = 1:iterations) %dopar% {
    
    boot_df <- draw_bootstrap_samples(dat)
    
    model <- mgcv::gam(formula = p_scaled ~ s(LEN_MIDPOINT, bs = "cr", k = k) + s(MATCHUP, bs = "re", by = dummy_var),
                       family = binomial(link = "logit"),
                       data = boot_df)
    
    fit_df <- data.frame(LEN_MIDPOINT = lengths,
                         MATCHUP = dat$MATCHUP[1],
                         dummy_var = 0) # random effects off
    
    fit_df$fit <- predict(model, newdata = fit_df, type = "link")
    
    return(fit_df)
    
  }
  
  doParallel::stopImplicitCluster()
  
  return(bootstrap_output)
  
}
