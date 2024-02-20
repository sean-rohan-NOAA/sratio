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

size = 10:50

abundance = round(dnorm(size, mean = 35, sd = 10) * 1e5 * (1-rnorm(length(size), mean = 0, sd = 0.1)))

selectivity_opts1 = list(type = "asymptotic",
                          begin_top = 20,
                          ln_sd1 = 5)

selectivity_opts2 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10)

# Setup selectivity
s_at_size1 <- selectivity_at_size(size = size, 
                                 selectivity_opts = selectivity_opts1)

s_at_size2 <- selectivity_at_size(size = size, 
                                 selectivity_opts = selectivity_opts2)

c_at_size1 <- catch_at_size(size = size, 
                            abundance = abundance, 
                            selectivity = s_at_size1,
                            n_size_samples = 10)

plot(size, s_at_size, type = 'l')
lines(size, s_at_size2, col = "red")


# Import numbers-at-length/age

# Setup proportion available at-length (numbers * normalized distribution)

# Draw from Dirichelet multinomial

# Conduct selectivity ratio analysis

# Conduct SCCAL analysis

# Save results