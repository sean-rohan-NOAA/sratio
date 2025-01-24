library(sratio)
library(selfisher)

# Densities: 
# ---- Equal: 5, 10, 20, 50; 
# Effort: 1 versus 0.5, 1 versus 1
# sratio: 2, 1.5, 1.1, 1, 0.9
# Selectivity l v l, dn vs l
# Pairs: 20, 50, 100

# 20 size classes
sim_settings <- 
  expand.grid(
    density = c(6, 10, 20, 50),
    effort1 = 1,
    effort2 = c(1, 0.5),
    q_ratio = c(2, 1.5, 1.1, 1),
    n_pairs = c(20, 50, 100),
    s_pattern1 = "logistic1",
    s_pattern2 = c("logistic1", "logistic2"),
    n_sizes = 20,
    n_bootstrap = 100,
    sratio_type = "absolute",
    # weighting_method = c("none", "count", "residuals_by_count"),
    stringsAsFactors = FALSE
  ) |>
  dplyr::filter(s_pattern1 == s_pattern2, 
                effort1 > effort2, 
                q_ratio == 1)

# Set maximum selectivity
max_s <- 0.9

# Flag scenarios where selectivity is identical
sim_settings <- sim_settings |>
  dplyr::filter(s_pattern1 == s_pattern2, q_ratio == q_ratio) |>
  dplyr::mutate(s_identical = TRUE) |>
  dplyr::full_join(sim_settings)

sim_settings$s_identical[is.na(sim_settings$s_identical)] <- FALSE

sim_settings$scenario <- 1:nrow(sim_settings)

# Create output list to store bootstrap results
bootstrap_results <- vector(mode = "list", 
                            length = nrow(sim_settings))

dir.create(here::here("analysis", "eval_sratio", "output", "bootstrap_samples_fits"),
           showWarnings = FALSE)


start_time <- Sys.time()

for(ii in 1:nrow(sim_settings)) {
  
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
  # n_draws <- sim_settings$n_bootstrap[ii]
  
  # Set selectivity at size as the product of selectivity and catchability ratio (q_ratio)
  sizes <- 1:sim_settings$n_sizes[ii]
  
 library(sratio)

# Densities: 
# ---- Equal: 5, 10, 20, 50; 
# Effort: 1 versus 0.5, 1 versus 1
# sratio: 2, 1.5, 1.1, 1, 0.9
# Selectivity l v l, dn vs l
# Pairs: 20, 50, 100

# 20 size classes
sim_settings <- 
  expand.grid(
    density = c(6, 10, 20, 50),
    effort1 = 1,
    effort2 = c(1, 0.5),
    q_ratio = c(2, 1.5, 1.1, 1),
    n_pairs = c(20, 50, 100),
    s_pattern1 = "logistic1",
    s_pattern2 = c("logistic1", "logistic2"),
    n_sizes = 20,
    n_bootstrap = 100,
    sratio_type = "absolute",
    weighting_method = NA,
    stringsAsFactors = FALSE
  ) |>
  dplyr::filter(s_pattern1 == s_pattern2, 
                effort1 > effort2, 
                q_ratio == 1)

# Set maximum selectivity
max_s <- 0.9

# Flag scenarios where selectivity is identical
sim_settings <- sim_settings |>
  dplyr::filter(s_pattern1 == s_pattern2, q_ratio == q_ratio) |>
  dplyr::mutate(s_identical = TRUE) |>
  dplyr::full_join(sim_settings)

sim_settings$s_identical[is.na(sim_settings$s_identical)] <- FALSE

sim_settings$scenario <- 1:nrow(sim_settings)

# Create output list to store bootstrap results
bootstrap_results <- vector(mode = "list", 
                            length = nrow(sim_settings))

dir.create(here::here("analysis", "eval_sratio", "output", "bootstrap_samples_fits"),
           showWarnings = FALSE)


start_time <- Sys.time()

for(ii in 1:nrow(sim_settings)) {

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
  
  set.seed(999)
  
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
  
  sf_mod <- selfisher(p12 ~ offset(log(qratio)) + bs(size, df = 4) + (1 | block), 
                      data = samples, 
                      total = total_count, 
                      haul = block, 
                      psplit = FALSE)
  
  newdata <- data.frame(size = 1:20,
                         total_count = 1,
                         block = NA,
                         qratio = samples$qratio[1]
  )
  
  predict(sf_mod, newdata = newdata)
  
  sf_boot <- selfisher::bootSel(sf_mod, 
                                nsim = 100,
                                type = "nonparametric", 
                                # parallel = "multicore",
                                # cl = 4,
                                FUN=function(sf_mod){predict(sf_mod, newdata=newdata, type = "ratio")})
  
  sf_quantiles <- apply(sf_boot$t, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
  newdata[,c("q025", "q250", "q500", "q750", "q975")] <- t(sf_quantiles)
  
  ggplot() +
    geom_path(data = size_selectivity_density,
              mapping = aes(x = size, y = s12), color = "red") +
    geom_ribbon(data = newdata,
                mapping = aes(x = size, ymin = q025, ymax = q975),
                alpha = 0.3) +
    geom_path(data = newdata,
              mapping = aes(x = size, y = q250), linetype = 3) +
    geom_path(data = newdata,
              mapping = aes(x = size, y = q500)) +
    geom_path(data = newdata,
              mapping = aes(x = size, y = q750), linetype = 3)

  
  error_summary <-  size_selectivity_density |>
    dplyr::select(size, s12_actual = s12) |>
    dplyr::inner_join(boot_fit, by = "size") |>
    dplyr::group_by(size) |>
    dplyr::summarise(s_12_mean_diff = mean(s12 - s12_actual),
                     s_12_median_diff = median(s12 - s12_actual),
                     s_12_mre = mean((s12-s12_actual)/s12_actual))
  
  # Calculate selectivity ratio quantiles
  bootstrap_quantiles <- boot_fit |>
    dplyr::group_by(size) |>
    dplyr::summarise(
      sratio_q025 = quantile(s12, 0.025),
      sratio_q250 = quantile(s12, 0.25),
      sratio_q500 = quantile(s12, 0.5),
      sratio_q750 = quantile(s12, 0.75),
      sratio_q975 = quantile(s12, 0.975)
    ) |>
    dplyr::mutate(type = "Bootstrap")
  
  bootstrap_results[[ii]] <- 
    size_selectivity_density |>
    dplyr::inner_join(bootstrap_quantiles, by = "size") |>
    dplyr::inner_join(error_summary, by = "size")
  
  end_time <- Sys.time()
  
  cat( "Time elapsed (", ii, "/", nrow(sim_settings), "): ", 
       round(difftime(end_time, start_time, units = "mins"), 2), 
       " minutes \n")
  
}

dir.create(here::here("analysis", "eval_sratio", "output"), recursive = TRUE, showWarnings = FALSE)

saveRDS(object = bootstrap_results,
        file = here::here("analysis", "eval_sratio", "output", "sim_bootstrap.rds"))

ggplot() +
  geom_ribbon(data = bootstrap_results[[ii]],
              mapping = aes(x = size,
                            ymin = sratio_q025,
                            max = sratio_q975),
              alpha = 0.5,
              fill = "grey20") +
  geom_path(data = bootstrap_results[[ii]],
            mapping = aes(x = size,
                          y = sratio_q250),
            linetype = 3) +
  geom_path(data = bootstrap_results[[ii]],
            mapping = aes(x = size,
                          y = sratio_q750),
            linetype = 3) +
  geom_path(data = bootstrap_results[[ii]],
            mapping = aes(x = size,
                          y = sratio_q500)) +
  geom_path(data = size_selectivity_density,
            mapping = aes(x = size, y = s12, color = "True value")) +
  scale_y_log10(name = expression(italic(S['L,30,15'])~(SR)), expand = c(0.05, 0.05)) +
  theme_bw()


bootstrap_results <- readRDS(file = here::here("analysis", "eval_sratio", "output", "sim_bootstrap.rds"))


bootstrap_df <- do.call("rbind", bootstrap_results)


# Example scenario

ex_scenario <- bootstrap_df |>
  dplyr::filter(s_pattern1 == s_pattern2, 
                effort1 > effort2, 
                q_ratio == 1) |>
  dplyr::filter(density == 10, n_pairs == 50, weighting_method == "residuals_by_count")



ggplot() +
  geom_ribbon(data = ex_scenario,
              mapping = aes(x = size, ymin = sratio_q025, ymax = sratio_q975, fill = "95% CI"),
              alpha = 0.5) +
  geom_path(data = ex_scenario,
            mapping = aes(x = size, y = sratio_q500, linetype = "Median"), 
            linewidth = 1.15) +
  geom_path(data = ex_scenario |>
              dplyr::select(size, sratio_q750, sratio_q250) |>
              tidyr::pivot_longer(cols = c("sratio_q750", "sratio_q250")),
            mapping = aes(x = size, y = value, linetype = "50% CI", group = name), 
            linewidth = 1.15) +
  geom_path(data = ex_scenario,
            mapping = aes(x = size, y = s12, color = "Actual"), 
            linewidth = 1.15) +
  ggtitle("Density = 10, N_hauls = 50, Weighting method = 'Residual by count (inverse variance)'") +
  scale_linetype_manual(values = c("Median" = 1, "50% CI" = 3)) +
  scale_fill_manual(values = c("95% CI" = "grey50")) +
  scale_color_manual(values = c("50% CI" = "black", "Median" = "black", "Actual" = "red")) +
  scale_x_continuous(name = "Size") +
  scale_y_continuous(name = expression('Selectivity ratio, '~S['1,2'])) +
  theme_bw() +
  theme(legend.title = element_blank())


# Mean relative error ----
max_mre <- 200

ggplot() +
  geom_path(data = dplyr::filter(bootstrap_df, 
                                 s_pattern1 == s_pattern2, 
                                 effort1 > effort2, 
                                 q_ratio == 1) |>
              dplyr::mutate(s_12_mre_pct = s_12_mre*100) |>
              dplyr::mutate(s_12_mre_pct = if_else(s_12_mre_pct > max_mre, max_mre, s_12_mre_pct)),
            mapping = aes(x = size, y = s_12_mre_pct, color = weighting_method, group = scenario),
            linewidth = 1.5) + 
  scale_y_continuous(name = "Mean relative error (%)", limits = c(-10, max_mre+1), expand = c(0,0)) +
  scale_x_continuous(limits = c(1,20)) +
  scale_color_colorblind(name = "Weighting method") +
  facet_grid(density ~ n_pairs) +
  ggtitle("Mean relative error (%) in selectivity ratio for hauls with equal catchability and selectivity (logistic), but where gear #2 has half the effort.\nThree different observation weighting methods.\nColumns = Number of hauls, Rows = Fish density.") +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot() +
  geom_path(data = dplyr::filter(bootstrap_df, 
                                 s_pattern1 == s_pattern2, 
                                 effort1 > effort2, 
                                 q_ratio == 1) |>
              dplyr::mutate(s_12_mre_pct = s_12_mre*100) |>
              dplyr::mutate(s_12_mre_pct = if_else(s_12_mre_pct > max_mre, max_mre, s_12_mre_pct)),
            mapping = aes(x = size, y = s_12_mre_pct, color = weighting_method, group = scenario),
            linewidth = 1.2) + 
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(name = "Mean relative error (%)", limits = c(-10, 10), expand = c(0,0)) +
  scale_x_continuous(limits = c(1,20)) +
  scale_color_colorblind(name = "Weighting method") +
  facet_grid(density ~ n_pairs) +
  ggtitle("Mean relative error (%) in selectivity ratio for hauls with equal catchability and selectivity (logistic), but where gear #2 has half the effort (for focal range of -10% to 10%.\nThree observation weighting methods.\nColumns = Number of hauls, Rows = Fish density.") +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot() +
  geom_bar(data = dplyr::filter(bootstrap_df, 
                                 s_pattern1 == s_pattern2, 
                                 effort1 > effort2, 
                                 q_ratio == 1) |>
              dplyr::select(size, n_pairs, density, mean_count1, mean_count2) |>
              unique() |>
              tidyr::pivot_longer(cols = c("mean_count1", "mean_count2")),
            mapping = aes(x = size, 
                          y = value, 
                          fill = name),
           position = "dodge",
           stat = "identity",
            linewidth = 1.2) + 
  scale_y_continuous(name = "Mean catch frequency", expand = c(0,0)) +
  scale_x_continuous(limits = c(1,20)) +
  scale_fill_colorblind(name = "Gear") +
  facet_grid(density ~ n_pairs) +
  ggtitle("Mean number of lengths per haul for hauls with equal catchability and selectivity (logistic), but where gear #2 has half the effort.\nColumns = Number of hauls, Rows = Fish density.") +
  theme_bw() +
  theme(legend.position = "bottom")
  
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
  
  set.seed(999)
  
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
  
  
  # Two-stage bootstrap draw of samples, where the number of samples drawn is equal to n_bootstrap
  boot_samples <- 
    sratio::two_stage_bootstrap(
      count1 = samples$count1,
      count2 = samples$count2,
      size1 = samples$size,
      size2 = samples$size,
      block1 = samples$block,
      block2 = samples$block,
      treatment_name1 = 1,
      treatment_name2 = 2,
      n_draws = n_draws
    )
  
  # Fill effort and sampling factors because they aren't propagated through two_stage_bootstrap (yet)
  for(kk in 1:length(boot_samples)) {
    
    boot_samples[[kk]]$effort <- sim_settings$effort1[ii]
    boot_samples[[kk]]$effort[boot_samples[[kk]]$treatment == 2] <- sim_settings$effort2[ii]
    boot_samples[[kk]]$sampling_factor <- 1
    
  }
  
  # Fit models to bootstrap samples
  boot_fit <- 
    sratio::sratio_fit_bootstrap(
      x = boot_samples,
      treatment_order = c(1,2),
      treatment_col = "treatment",
      count_col = "count",
      size_col = "size",
      block_col = "new_block",
      effort_col = "effort",
      sampling_factor_col = "sampling_factor",
      gam_family = "binomial",
      gam_formula = p ~ s(size, bs = "tp", k = 8) + s(block, bs = "re") + offset(log_qratio),
      k = 8,
      obs_weight_control = obs_weight_control,
      sratio_type = sratio_type,
      n_cores = 4
    )
  
  saveRDS(object = list(boot_samples = boot_samples,
                        boot_fit = boot_fit),
          file = here::here("analysis", "eval_sratio", "output", "bootstrap_samples_fits",  paste0("bootstrap_sample_", ii, ".rds")))
  
  error_summary <-  size_selectivity_density |>
    dplyr::select(size, s12_actual = s12) |>
    dplyr::inner_join(boot_fit, by = "size") |>
    dplyr::group_by(size) |>
    dplyr::summarise(s_12_mean_diff = mean(s12 - s12_actual),
                     s_12_median_diff = median(s12 - s12_actual),
                     s_12_mre = mean((s12-s12_actual)/s12_actual))
  
  # Calculate selectivity ratio quantiles
  bootstrap_quantiles <- boot_fit |>
    dplyr::group_by(size) |>
    dplyr::summarise(
      sratio_q025 = quantile(s12, 0.025),
      sratio_q250 = quantile(s12, 0.25),
      sratio_q500 = quantile(s12, 0.5),
      sratio_q750 = quantile(s12, 0.75),
      sratio_q975 = quantile(s12, 0.975)
    ) |>
    dplyr::mutate(type = "Bootstrap")
  
  bootstrap_results[[ii]] <- 
    size_selectivity_density |>
    dplyr::inner_join(bootstrap_quantiles, by = "size") |>
    dplyr::inner_join(error_summary, by = "size")
  
  end_time <- Sys.time()
  
  cat( "Time elapsed (", ii, "/", nrow(sim_settings), "): ", 
       round(difftime(end_time, start_time, units = "mins"), 2), 
       " minutes \n")
  
}

dir.create(here::here("analysis", "eval_sratio", "output"), recursive = TRUE, showWarnings = FALSE)

saveRDS(object = bootstrap_results,
        file = here::here("analysis", "eval_sratio", "output", "sim_bootstrap.rds"))

ggplot() +
  geom_ribbon(data = bootstrap_results[[ii]],
              mapping = aes(x = size,
                            ymin = sratio_q025,
                            max = sratio_q975),
              alpha = 0.5,
              fill = "grey20") +
  geom_path(data = bootstrap_results[[ii]],
            mapping = aes(x = size,
                          y = sratio_q250),
            linetype = 3) +
  geom_path(data = bootstrap_results[[ii]],
            mapping = aes(x = size,
                          y = sratio_q750),
            linetype = 3) +
  geom_path(data = bootstrap_results[[ii]],
            mapping = aes(x = size,
                          y = sratio_q500)) +
  geom_path(data = size_selectivity_density,
            mapping = aes(x = size, y = s12, color = "True value")) +
  scale_y_log10(name = expression(italic(S['L,30,15'])~(SR)), expand = c(0.05, 0.05)) +
  theme_bw()


bootstrap_results <- readRDS(file = here::here("analysis", "eval_sratio", "output", "sim_bootstrap.rds"))


bootstrap_df <- do.call("rbind", bootstrap_results)


# Example scenario

ex_scenario <- bootstrap_df |>
  dplyr::filter(s_pattern1 == s_pattern2, 
                effort1 > effort2, 
                q_ratio == 1) |>
  dplyr::filter(density == 10, n_pairs == 50, weighting_method == "residuals_by_count")



ggplot() +
  geom_ribbon(data = ex_scenario,
              mapping = aes(x = size, ymin = sratio_q025, ymax = sratio_q975, fill = "95% CI"),
              alpha = 0.5) +
  geom_path(data = ex_scenario,
            mapping = aes(x = size, y = sratio_q500, linetype = "Median"), 
            linewidth = 1.15) +
  geom_path(data = ex_scenario |>
              dplyr::select(size, sratio_q750, sratio_q250) |>
              tidyr::pivot_longer(cols = c("sratio_q750", "sratio_q250")),
            mapping = aes(x = size, y = value, linetype = "50% CI", group = name), 
            linewidth = 1.15) +
  geom_path(data = ex_scenario,
            mapping = aes(x = size, y = s12, color = "Actual"), 
            linewidth = 1.15) +
  ggtitle("Density = 10, N_hauls = 50, Weighting method = 'Residual by count (inverse variance)'") +
  scale_linetype_manual(values = c("Median" = 1, "50% CI" = 3)) +
  scale_fill_manual(values = c("95% CI" = "grey50")) +
  scale_color_manual(values = c("50% CI" = "black", "Median" = "black", "Actual" = "red")) +
  scale_x_continuous(name = "Size") +
  scale_y_continuous(name = expression('Selectivity ratio, '~S['1,2'])) +
  theme_bw() +
  theme(legend.title = element_blank())


# Mean relative error ----
max_mre <- 200

ggplot() +
  geom_path(data = dplyr::filter(bootstrap_df, 
                                 s_pattern1 == s_pattern2, 
                                 effort1 > effort2, 
                                 q_ratio == 1) |>
              dplyr::mutate(s_12_mre_pct = s_12_mre*100) |>
              dplyr::mutate(s_12_mre_pct = if_else(s_12_mre_pct > max_mre, max_mre, s_12_mre_pct)),
            mapping = aes(x = size, y = s_12_mre_pct, color = weighting_method, group = scenario),
            linewidth = 1.5) + 
  scale_y_continuous(name = "Mean relative error (%)", limits = c(-10, max_mre+1), expand = c(0,0)) +
  scale_x_continuous(limits = c(1,20)) +
  scale_color_colorblind(name = "Weighting method") +
  facet_grid(density ~ n_pairs) +
  ggtitle("Mean relative error (%) in selectivity ratio for hauls with equal catchability and selectivity (logistic), but where gear #2 has half the effort.\nThree different observation weighting methods.\nColumns = Number of hauls, Rows = Fish density.") +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot() +
  geom_path(data = dplyr::filter(bootstrap_df, 
                                 s_pattern1 == s_pattern2, 
                                 effort1 > effort2, 
                                 q_ratio == 1) |>
              dplyr::mutate(s_12_mre_pct = s_12_mre*100) |>
              dplyr::mutate(s_12_mre_pct = if_else(s_12_mre_pct > max_mre, max_mre, s_12_mre_pct)),
            mapping = aes(x = size, y = s_12_mre_pct, color = weighting_method, group = scenario),
            linewidth = 1.2) + 
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(name = "Mean relative error (%)", limits = c(-10, 10), expand = c(0,0)) +
  scale_x_continuous(limits = c(1,20)) +
  scale_color_colorblind(name = "Weighting method") +
  facet_grid(density ~ n_pairs) +
  ggtitle("Mean relative error (%) in selectivity ratio for hauls with equal catchability and selectivity (logistic), but where gear #2 has half the effort (for focal range of -10% to 10%.\nThree observation weighting methods.\nColumns = Number of hauls, Rows = Fish density.") +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot() +
  geom_bar(data = dplyr::filter(bootstrap_df, 
                                s_pattern1 == s_pattern2, 
                                effort1 > effort2, 
                                q_ratio == 1) |>
             dplyr::select(size, n_pairs, density, mean_count1, mean_count2) |>
             unique() |>
             tidyr::pivot_longer(cols = c("mean_count1", "mean_count2")),
           mapping = aes(x = size, 
                         y = value, 
                         fill = name),
           position = "dodge",
           stat = "identity",
           linewidth = 1.2) + 
  scale_y_continuous(name = "Mean catch frequency", expand = c(0,0)) +
  scale_x_continuous(limits = c(1,20)) +
  scale_fill_colorblind(name = "Gear") +
  facet_grid(density ~ n_pairs) +
  ggtitle("Mean number of lengths per haul for hauls with equal catchability and selectivity (logistic), but where gear #2 has half the effort.\nColumns = Number of hauls, Rows = Fish density.") +
  theme_bw() +
  theme(legend.position = "bottom")
