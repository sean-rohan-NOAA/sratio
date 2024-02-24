#' Generate a normalized vector from normal, skew-normal or bimodal distributions
#' 
#' @param size Numeric vector of sizes
#' @param demographic_comp_pars A list containing values for the distribution
#' @export

scaled_distribution <- function(size, demographic_comp_pars = NULL) {
  
  stopifnot("scaled_distribution: demographic_comp_pars must be a list containing an object named 'distribution'" = "distribution" %in% names(demographic_comp_pars))
  
  if(demographic_comp_pars$distribution == "none") {
    
    dens <- rep(1, length(size))
    
  }
  
  if(demographic_comp_pars$distribution == "normal") {
    
    stopifnot("scaled_distribution: When demographic_comp_pars$distribution is normal, demographic_comp_pars must contain 1L numerical vectors of distribution parameters named 'mean' and 'sd'" = all(c("mean", "sd") %in% names(demographic_comp_pars)))
    
    dens <- sn::dsn(x = size,
                    xi = demographic_comp_pars$mean,
                    omega = demographic_comp_pars$sd,
                    alpha = 0)
    
  }
  
  if(demographic_comp_pars$distribution == "sn") {
    
    stopifnot("scaled_distribution: When demographic_comp_pars$distribution is 'sn', demographic_comp_pars must contain 1L numerical vectors of distribution parameters named 'xi', 'omega', and 'alpha'" = all(c("xi", "omega", "alpha") %in% names(demographic_comp_pars)))
    
    dens <- sn::dsn(x = size,
                    xi = demographic_comp_pars$xi,
                    omega = demographic_comp_pars$omega,
                    alpha = demographic_comp_pars$alpha)
    
  }
  
  if(demographic_comp_pars$distribution == "bimodal") {
    
    stopifnot("scaled_distribution: When demographic_comp_pars$distribution is 'bimodal', demographic_comp_pars must contain 1L numerical vectors of distribution parameters named 'mean1', 'sd1', 'mean2', and 'sd2'" = all(c("mean1", "sd1", "mean2", "sd2") %in% names(demographic_comp_pars)))
    
    dens <- sn::dsn(x = size,
                    xi = demographic_comp_pars$mean1,
                    omega = demographic_comp_pars$sd1,
                    alpha = 0) +
      sn::dsn(x = size,
              xi = demographic_comp_pars$mean2,
              omega = demographic_comp_pars$sd2,
              alpha = 0)
    
  }
  
  prop <- dens / max(dens)
  
  return(prop)
  
}