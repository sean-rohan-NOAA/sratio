#' Simulate catch for paired sample between two gears/methods
#' 
#' Simulate catch-at-size for paired sample between two gears/methods based on numbers-at-size, availability, local demographic structure, gear efficiency, effort (e.g. area swept), and selectivity.
#' 
#' @param size Sizes as a numeric vector
#' @param abundance Numbers-at-size as a numeric vector
#' @param availability Proportion of the total population that is within the sampled area.
#' @param demographic_comp_pars Local demographic composition parameters that are used to scale abundance as a list. Required named elements of the list vary depending on the distributional assumption (see ?scaled_distribution.).
#' @param gear_q1 Gear efficiency for gear/method #1 as a 1L numeric vector
#' @param gear_q2 Same as above for gear/method #2.
#' @param effort1 Effort for gear/method #1 as a 1L numeric vector
#' @param effort2 Same as above for gear/method #2.
#' @param selectivity_opts1 Selectivity model parameters for gear/method #1 as a list. Required named elements of the list vary depending on the selectivity assumption (see ?selectivity_at_size)
#' @param selectivity_opts2 Same as above for gear/method #2.
#' @param return_vars If TRUE, returns a list containing variable values used for the simulation. Otherwise, only returns the samples.
#' @returns A list containing size samples for both gears.
#' @examples
#' # A vector of sizes
#' size <- 10:55
#' 
#' # A vector of abundance-at-size
#' abundance = round(dnorm(size, mean = 35, sd = 10) * 1e5 * (1-rnorm(length(size), mean = 0, sd = 0.1)))
#' 
#' # Proportion of the total population that is available to a haul ----
#' availability = 0.001
#' 
#' demographic_comp_pars = list(distribution = "normal", mean = 30, sd = 10)
#' 
#' # Gear efficiency ----
#' gear_q1 = 1
#' gear_q2 = 0.8
#' 
#' # Effort for the treatments ----
#' effort1 = 0.5
#' effort2 = 1
#' 
#' # Options for the selectivity function ----
#' selectivity_opts1 = list(type = "asymptotic",
#'                          begin_top = 35,
#'                          ln_sd1 = 10)
#' 
#' selectivity_opts2 = list(type = "asymptotic",
#'                          begin_top = 35,
#'                          ln_sd1 = 10)
#' 
#' sample1 <- simulate_paired_sample(size = size, 
#'                                   abundance = abundance, 
#'                                   availability = availability, 
#'                                   demographic_comp_pars = demographic_comp_pars, 
#'                                   gear_q1 = gear_q1, 
#'                                   gear_q2 = gear_q2, 
#'                                   effort1 = effort1, 
#'                                   effort2 = effort2, 
#'                                   selectivity_opts1 = selectivity_opts1, 
#'                                   selectivity_opts2 = selectivity_opts2)
#' @export

simulate_paired_sample <- function(size, abundance, availability, demographic_comp_pars, gear_q1, gear_q2, effort1, effort2, selectivity_opts1, selectivity_opts2, return_vars = TRUE) {
  
  # Setup local abundance based on the proportional composition of the population and availability ----
  demographic_proportion <- scaled_distribution(size = size, 
                                                demographic_comp_pars = demographic_comp_pars)
  
  demographic_proportion <- demographic_proportion / max(demographic_proportion)
  
  if(class(availability) == "list") {
    availability <- draw_from_dist(x = availability)
  }
  
  local_density <- abundance * demographic_proportion * availability
  
  # Setup selectivity for two gears -----
  s_at_size1 <- selectivity_at_size(size = size, 
                                    selectivity_opts = selectivity_opts1)
  
  s_at_size2 <- selectivity_at_size(size = size, 
                                    selectivity_opts = selectivity_opts2)
  
  # Total catch in numbers with rounding of fractional values (scales with density) ----
  sample_size1 <- sum(local_density * s_at_size1 * effort1 * gear_q1)
  
  sample_size1 <- sample(x = c(0,1), 
                         size = 1, 
                         prob = c(1-(sample_size1 %% 1), sample_size1 %% 1)) + floor(sample_size1)
  
  sample_size2 <- sum(local_density * s_at_size2 * effort2 * gear_q2)
  
  sample_size2 <- sample(x = c(0,1), 
                         size = 1, 
                         prob = c(1-(sample_size2 %% 1), sample_size2 %% 1)) + floor(sample_size2)
  
  # Catch-at-size based on Dirichlet-multinomial draws from the product of local abundance and selectivity ----
  c_at_size1 <- catch_at_size(size = size, 
                              abundance = local_density, 
                              selectivity = s_at_size1,
                              n_size_samples = sample_size1)
  
  c_at_size2 <- catch_at_size(size = size, 
                              abundance = local_density, 
                              selectivity = s_at_size2,
                              n_size_samples = sample_size2)
  
  c_at_size1$effort <- effort1
  c_at_size2$effort <- effort2
  c_at_size1$treatment <- 1
  c_at_size2$treatment <- 2
  
  samples <- rbind(c_at_size1, c_at_size2)
  
  if(return_vars) {
    output <- list(samples = samples, 
                   vars = list(availability = availability,
                       demographic_comp_pars = demographic_proportion,
                       local_density = local_density,
                       gear_q1 = gear_q1,
                       gear_q2 = gear_q1, 
                       effort1 = effort1, 
                       effort2 = effort2,
                       selectivity1 = s_at_size1, 
                       selectivity2 = s_at_size2))
  } else {
    
    output <- list(samples = samples)
    
  }
  
  return(output)
  
}