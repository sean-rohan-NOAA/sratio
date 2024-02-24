# Selectivity simulation

# Parameters
# - Sample size (hauls)
# - Selectivity function
# - Size vector
# - Numbers-at-size
# - Demographic distribution
# - Sample sizes/proportion
# - Size bins for analysis
# - Bootstrap draws per iteration
# - Number of iterations
# - Sample size (hauls)

# Outputs
# - Bootstrap samples (temporary)
# - Draw index
# - Settings
# - Fits (by draw)
library(sratio)

# Function arguments ----

demographic_comp_pars = list(distribution = "normal", mean = 30, sd = 10)
# demographic_comp_pars <- list(distribution = "sn", xi = 30, omega = 10, alpha = 3)


# A vector of sizes
size <- 10:55

# A vector of abundance-at-size
abundance = round(dnorm(size, mean = 35, sd = 10) * 1e5 * (1-rnorm(length(size), mean = 0, sd = 0.1)))

# Proportion of the total population that is available to a haul ----
# availability = 0.0003

availability = list(distribution = "normal", mean = 0.0004, sd = 0.00002)

demographic_comp_pars = list(distribution = "normal", mean = 30, sd = 10)

# Gear efficiency ----
gear_q1 = 1
gear_q2 = 0.8

# Effort for the treatments ----
effort1 = 0.5
effort2 = 1

selectivity_opts1 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10) 

selectivity_opts2 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10)

n_pairs <- 40


for(ii in 1:n_pairs) {
  
}
  
  sample1 <- simulate_paired_sample(size = size, 
                                    abundance = abundance, 
                                    availability = availability, 
                                    demographic_comp_pars = demographic_comp_pars, 
                                    gear_q1 = gear_q1, 
                                    gear_q2 = gear_q2, 
                                    effort1 = effort1, 
                                    effort2 = effort2, 
                                    selectivity_opts1 = selectivity_opts1, 
                                    selectivity_opts2 = selectivity_opts2)

  
}

# Options for the selectivity function ----
selectivity_opts1 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10)

selectivity_opts2 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10)




plot(size, s_at_size1*gear_q1, type = 'l')
lines(size, s_at_size2*gear_q2, col = "red")


# Draw from Dirichelet multinomial

# Conduct selectivity ratio analysis

# Conduct SCCAL analysis

# Save results