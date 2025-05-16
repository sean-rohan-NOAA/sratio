library(sratio)
library(gearcalib)

# Densities: 
# ---- Equal: 6, 10, 20, 50; 
# Effort: 1 versus 0.5, 1 versus 1
# sratio: 2, 1.5, 1.1, 1, 0.9
# Selectivity l v l, dn vs l
# Pairs: 20, 50, 100

ii = 1

# 20 size classes
sim_settings <- 
  expand.grid(
    density = c(6, 10, 20, 50, 100),
    effort1 = 1,
    effort2 = 0.5,
    q_ratio = 1,
    n_pairs = c(20, 50, 100),
    s_pattern1 = "logistic1",
    s_pattern2 = "logistic1",
    n_sizes = 13,
    n_bootstrap = 100,
    sratio_type = "absolute",
    weighting_method = c("none", "count", "residuals_by_count"),
    stringsAsFactors = FALSE
  ) |>
  dplyr::filter(s_pattern1 == s_pattern2, 
                effort1 > effort2, 
                q_ratio == 1)

# Set maximum selectivity
max_s <- 0.9

sim_settings$scenario <- 1:nrow(sim_settings)

# Create output list to store bootstrap results
bootstrap_results <- vector(mode = "list", 
                            length = nrow(sim_settings))

thygesen_results <- vector(mode = "list", 
                           length = nrow(sim_settings))

dir.create(here::here("analysis", "eval_sratio", "output", "bootstrap_samples_fits"),
           showWarnings = FALSE)


start_time <- Sys.time()

# Type of selectivity ratio to calculate: absolute or relative
sratio_type <-  sim_settings$sratio_type[ii]

# Set observation weighting method for models
obs_weight_control <- 
  switch(
    sim_settings$weighting_method[ii],
    "none" = 
      list(method = "none",
           max_count = Inf,
           residual_type = NA,
           normalize_weights = FALSE),
    "count" = 
      list(method = "count", 
           max_count = Inf,
           residual_type = NA,
           normalize_weights = FALSE),
    "residuals_by_count" = 
      list(method = "residuals_by_count", 
           max_count = Inf,
           residual_type = "absolute",
           normalize_weights = FALSE)
  )

# Number of bootstrap samples to draw to estimate CIs
n_draws <- sim_settings$n_bootstrap[ii]

# Set selectivity at size as the product of selectivity and catchability ratio (q_ratio)
sizes <- 1:sim_settings$n_sizes[ii]

s1_opts <- switch(
  sim_settings$s_pattern1[ii],
  "logistic1" = list(type = "logistic",
                     begin_top = 7,
                     ln_sd1 = 6),
  "logistic2" = list(type = "logistic",
                     begin_top = 10,
                     ln_sd1 = 6),
  "doublenormal1" = doublenormal1
)

s2_opts <- switch(
  sim_settings$s_pattern2[ii],
  "logistic1" = list(type = "logistic",
                     begin_top = 7,
                     ln_sd1 = 6),
  "logistic2" = list(type = "logistic",
                     begin_top = 10,
                     ln_sd1 = 6),
  "doublenormal1" = doublenormal1
)

s1 <- selectivity_at_size(
  size = sizes, 
  selectivity_opts = s1_opts
)

s2 <- selectivity_at_size(
  size = sizes, 
  selectivity_opts = s2_opts
)

# Catchability ratio adjustment
s2 <- s2 / sim_settings$q_ratio[ii]

# Set maximum selectivity
s2[s2 > max_s] <- max_s
s1[s1 > max_s] <- max_s


s12 <- s1/s2

# Number of encounters for each haul is the product of effort and density
encounter1 <- sim_settings$density[ii] * sim_settings$effort1[ii]

encounter2 <- sim_settings$density[ii] * sim_settings$effort2[ii]

# Settings for the iith simulation
size_selectivity_density <- 
  data.frame(
    scenario = sim_settings$scenario[ii],
    size = sizes,
    density = sim_settings$density[ii],
    effort1 = sim_settings$effort1[ii],
    effort2 = sim_settings$effort2[ii],
    encounter1 = encounter1,
    encounter2 = encounter2,
    q_ratio = sim_settings$q_ratio[ii],
    s1 = s1,
    s2 = s2,
    s12 = s12,
    n_pairs = sim_settings$n_pairs[ii],
    s_pattern1 = sim_settings$s_pattern1[ii],
    s_pattern2 = sim_settings$s_pattern2[ii],
    weighting_method = sim_settings$weighting_method[ii]
  )


# Draw samples for each haul based on the number of encounters and retention probability at size
# Retention for each encounter is randomly drawn from a binomial distribution. 
# The retention probability for each encounter is the product of q and selectivity for the size class.
# The number of samples drawn is equal to the number of paired samples.

# It is assumed that all individuals caught are lengthed (sampling_factor1 = sampling_factor2 = 1)
samples <- data.frame()

for(jj in 1:nrow(size_selectivity_density)) {
  
  len_sample <- data.frame(
    length = size_selectivity_density$size[jj],
    n1 = rbinom(n = size_selectivity_density$n_pairs[jj],
                size = size_selectivity_density$encounter1[jj], 
                prob = size_selectivity_density$s1[jj]),
    n2 = rbinom(n = size_selectivity_density$n_pairs[jj], 
                size = size_selectivity_density$encounter2[jj], 
                prob = size_selectivity_density$s2[jj])
  )
  
  # Calculate the selectivity ratio for each haul
  len_sample <- sratio::selectivity_ratio(
    size = len_sample$length,
    count1 = len_sample$n1,
    count2 = len_sample$n2,
    effort1 = size_selectivity_density$effort1[jj],
    effort2 = size_selectivity_density$effort2[jj],
    sampling_factor1 = 1,
    sampling_factor2 = 1,
    sratio_type = sratio_type
  )
  
  # Setup blocks for random effects
  len_sample$block <- 1:nrow(len_sample)
  
  samples <- len_sample |>
    dplyr::bind_rows(samples)
  
}

# Fill total count (used for weighting methods)
samples$total_count <- samples$count1 + samples$count2

samples <- dplyr::filter(samples, !is.na(p12))

size_selectivity_density <- 
  dplyr::inner_join(
    size_selectivity_density,
    samples |>
      dplyr::group_by(size) |>
      dplyr::summarise(mean_count1 = mean(count1),
                       mean_count2 = mean(count2),
                       median_count1 = median(count1),
                       median_count2 = median(count2),
                       n_positive_catch1 = sum(count1 > 0),
                       n_positive_catch2 = sum(count2 > 0))
  )

###

test <-
  samples |>
  dplyr::select(size, block, count1, count2) |>
  tidyr::pivot_longer(cols = c("count1", "count2"), values_to = "count") |>
  dplyr::arrange(size, block) |>
  tidyr::pivot_wider(names_from = size, values_from = count, values_fill = 0)

test <-
  data.frame(name = c("count1", "count2"),
             effort = c(1, 0.5)) |>
  dplyr::inner_join(
    test
  )

d_input <-
  list(
    N = as.matrix(test[, 5:ncol(test)]),
    SweptArea = test$effort,
    group = factor(test$block),
    Gear = factor(test$name)
  )

out <- gearcalibFitNB(d = d_input, model = "poisson")

plot.gearcalibFit(out)
