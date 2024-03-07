# Selectivity simulation

library(sratio)

# Number of cores for parallel processing ----
n_cores = 4

# Sampling simulation parameters ----
# - Number of realizations
n_sims = 50

# - Sample size (hauls)
n_pairs = 40

# - Selectivity function
selectivity_opts1 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10) 

selectivity_opts2 = list(type = "asymptotic",
                         begin_top = 35,
                         ln_sd1 = 10)

# - Gear efficiency
gear_q1 = 1
gear_q2 = 0.8

# - Sampling effort (i.e., area swept)
effort1 = 0.5
effort2 = 1

# - Fish size vector
size = 10:55

# - Numbers-at-size
abundance = round(dnorm(size, mean = 35, sd = 10) * 1e5 * (1-rnorm(length(size), mean = 0, sd = 0.1)))

# - Local demographic composition
demographic_comp_pars = list(distribution = "normal", mean = 30, sd = 10)

# - Local availability
availability = list(distribution = "normal", mean = 0.0004, sd = 0.00008)

# Bootstrap parameters ----
# - Size bin width for models 
size_bin_width = 5

# - Number of bootstrap samples
n_boot_samples = 1000

# GAM knots
gam_knots = 4

# Minimum sample size
min_sample_size = 10

# Method for scaling catch comparison rate for beta regression (see ?sratio::scale_for_betareg)
scale_method <- "sc"


settings = list(n_sims = n_sims,
                n_pairs = n_pairs, 
                selectivity_opts1 = selectivity_opts1,
                selectivity_opts2 = selectivity_opts2,
                gear_q1 = gear_q1,
                gear_q2 = gear_q2,
                effort1 = effort1,
                effort2 = effort2,
                size = size,
                abundance = abundance,
                demographic_comp_pars = demographic_comp_pars,
                availability = availability,
                size_bin_width = size_bin_width,
                n_boot_samples = n_boot_samples,
                gam_knots = gam_knots,
                min_sample_size = min_sample_size,
                scale_method = scale_method,
                n_cores = n_cores)

run_selectivity_simulation <- function(settings) {
  
  list2env(settings, envir = environment())
  
  # Draw samples
  catch_samples <- vector(mode = "list", length = n_sims)
  
  cat("Drawing catch samples\n")
  for(ii in 1:n_sims) {
    
    paired_samples <- data.frame()
    
    for(jj in 1:n_pairs) {
      
      draw <- simulate_paired_sample(size = size, 
                                     abundance = abundance, 
                                     availability = availability, 
                                     demographic_comp_pars = demographic_comp_pars, 
                                     gear_q1 = gear_q1, 
                                     gear_q2 = gear_q2, 
                                     effort1 = effort1, 
                                     effort2 = effort2, 
                                     selectivity_opts1 = selectivity_opts1, 
                                     selectivity_opts2 = selectivity_opts2,
                                     return_vars = FALSE)$samples
      
      draw$haul <- jj
      
      paired_samples <- rbind(paired_samples, draw)
      
    }
    
    hauls_with_samples <- dplyr::group_by(paired_samples, haul, treatment) |>
      dplyr::summarise(total = n(), .groups = "keep") |>
      dplyr::filter(total >= min_sample_size) |>
      dplyr::select(-total) |>
      dplyr::group_by(haul) |>
      dplyr::summarise(total = n(), .groups = "keep") |>
      dplyr::filter(total == 2) |>
      dplyr::select(-total) |>
      as.data.frame()
    
    paired_samples <- dplyr::group_by(paired_samples, haul, treatment, size, effort) |>
      dplyr::summarise(count = n(), .groups = "keep") |>
      dplyr::filter(haul %in% hauls_with_samples$haul) |>
      as.data.frame()
    
    # Filter based on minimum sample size
    catch_samples[[ii]] <- paired_samples
    
  }
  
  
  # Draw bootstrap samples, fit models, and generate predictions
  boot_samples <- vector(mode = "list", length = n_sims)
  boot_results <- vector(mode = "list", length = n_sims)
  boot_fit_sratio <- vector(mode = "list", length = n_sims)
  
  start_time <- Sys.time()
  for(kk in 1:n_sims) {
    
    cat(paste0("Running analysis on sample ", kk, " out of ", n_sims, "\n"))
    
    sel_sim <- catch_samples[[kk]]
    
    sel_sim$size_bin <- round_any(sel_sim$size, accuracy = size_bin_width)
    
    boot_samples[[kk]] <- sratio::two_stage_bootstrap(count1 = sel_sim$n[sel_sim$treatment == 1],
                                                      count2 = sel_sim$n[sel_sim$treatment == 2],
                                                      size1 = sel_sim$size_bin[sel_sim$treatment == 1],
                                                      size2 = sel_sim$size_bin[sel_sim$treatment == 2],
                                                      block1 = sel_sim$haul[sel_sim$treatment == 1],
                                                      block2 = sel_sim$haul[sel_sim$treatment == 2],
                                                      n_draws = n_boot_samples,
                                                      seed = NULL,
                                                      treatment_name1 = 1,
                                                      treatment_name2 = 2)
    
    # Add effort to bootstrap samples
    for(ll in 1:length(boot_samples[[kk]])) {
      
      boot_samples[[kk]][[ll]] <- dplyr::inner_join(boot_samples[[kk]][[ll]],
                                                    data.frame(treatment = c(1, 2),
                                                               effort = c(effort1, effort2)),
                                                    by = "treatment")
      
    }
    
    # Selectivity ratio model selection between binomial and beta binomial based on matchup-level cross validationt
    samples_wide <- catch_samples[[kk]] |>
      tidyr::pivot_wider(names_from = "treatment", 
                         values_from = "count", 
                         values_fill = 0, 
                         names_prefix = "count") |>
      dplyr::mutate(effort1 = effort1,
                    effort2 = effort2)
    
    pratio <- sratio_cv(count1 = samples_wide$count1,
                        count2 = samples_wide$count2,
                        effort1 = samples_wide$effort1,
                        effort2 = samples_wide$effort2,
                        size = samples_wide$size,
                        block = samples_wide$haul,
                        k = gam_knots, # GAM knots
                        n_cores = n_cores,
                        scale_method = "sv")
    
    best_model <- data.frame(model = c("binomial", "beta"), 
                             rmse = c(sqrt(mean((pratio$cv_fit_logit-pratio$p)^2)),
                                      sqrt(mean((pratio$cv_fit_beta-pratio$p)^2))),
                             gam_knots = gam_knots) |>
      dplyr::mutate(best = rmse == min(rmse))
    
    best_model <- best_model$model[best_model$best][1]
    
    boot_sratio <- sratio::sratio_fit_bootstrap(x = boot_samples[[kk]],
                                                treatment_order = c(1, 2),
                                                size_col = "size",
                                                block_col = "block",
                                                treatment_col = "treatment",
                                                count_col = "count",
                                                effort_col = "effort",
                                                gam_family = best_model,
                                                k = gam_knots,
                                                scale_method = scale_method,
                                                n_cores = n_cores)
    
    boot_fit_sratio[[kk]] <- boot_fit_sratio
    
    boot_results[[kk]] <- boot_sratio |>
      dplyr::group_by(size) |>
      dplyr::summarise(sratio_q025 = quantile(s21, 0.025),
                       sratio_q250 = quantile(s21, 0.25),
                       sratio_q500 = quantile(s21, 0.5),
                       sratio_q750 = quantile(s21, 0.75),
                       sratio_q975 = quantile(s21, 0.975)) |>
      dplyr::mutate(true_s1 = selectivity_at_size(size = size, selectivity_opts = selectivity_opts1) * gear_q1,
                    true_s2 = selectivity_at_size(size = size, selectivity_opts = selectivity_opts2) * gear_q2,
                    true_sratio = true_s1/true_s2,
                    draw = kk)
    
  }
  
  output <- list(bootstrap_results = boot_results,
                 bootstrap_fits = boot_fit_sratio,
                 best_model = best_model,
                 catch_samples = catch_samples,
                 settings = list(n_sims = n_sims,
                                 n_pairs = n_pairs, 
                                 selectivity_opts1 = selectivity_opts1,
                                 selectivity_opts2 = selectivity_opts2,
                                 gear_q1 = gear_q1,
                                 gear_q2 = gear_q2,
                                 effort1 = effort1,
                                 effort2 = effort2,
                                 size = size,
                                 abundance = abundance,
                                 demographic_comp_pars = demographic_comp_pars,
                                 availability = availability,
                                 size_bin_width = size_bin_width,
                                 n_boot_samples = n_boot_samples,
                                 gam_knots = gam_knots,
                                 min_sample_size = min_sample_size,
                                 scale_method = scale_method,
                                 n_cores = n_cores)
  )
  
  end_time <- Sys.time()
  
  cat(paste0("Total run time: ", end_time - start_time))
  
  return(output)
  
}


test <- run_selectivity_simulation(settings = settings)


test <- do.call(what = "rbind", boot_results)

test_plot <- dplyr::mutate(test,
                           resid = true_sratio-sratio_q500)

ggplot() +
  geom_path(data = dplyr::group_by(test_plot, size) |>
              dplyr::summarise(mean_resid = mean(resid, na.rm = TRUE),
                               n = n()),
            mapping = aes(x = size, y = mean_resid))


# - True selectivity at size


ggplot() +
  geom_ribbon(data = bootstrap_quantiles,
              mapping = aes(x = size,
                            ymin = p_q025,
                            max = p_q975),
              alpha = 0.5,
              fill = "grey20") +
  geom_path(data = bootstrap_quantiles,
            mapping = aes(x = size,
                          y = p_q250),
            linetype = 3) +
  geom_path(data = bootstrap_quantiles,
            mapping = aes(x = size,
                          y = p_q750),
            linetype = 3) +
  geom_path(data = bootstrap_quantiles,
            mapping = aes(x = size,
                          y = p_q500)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_y_continuous(name = expression(italic(p['L,30,15'])), limits = c(0,1)) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw()



ggplot() +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(data = bootstrap_quantiles,
              mapping = aes(x = size,
                            ymin = sratio_q025,
                            max = sratio_q975),
              alpha = 0.5,
              fill = "grey20") +
  geom_path(data = bootstrap_quantiles,
            mapping = aes(x = size,
                          y = sratio_q250),
            linetype = 3) +
  geom_path(data = bootstrap_quantiles,
            mapping = aes(x = size,
                          y = sratio_q750),
            linetype = 3) +
  geom_path(data = bootstrap_quantiles,
            mapping = aes(x = size,
                          y = sratio_q500)) +
  geom_path(mapping = aes(x = size, y = sratio), color = "red") +
  scale_y_log10(name = expression(italic(S['L,15,30'])~(SR)), expand = c(0.05, 0.05)) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw()

