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
n_boot_samples = 100

# GAM knots
gam_knots = 5

# Minimum sample size
min_sample_size = 10

# Method for scaling catch comparison rate for beta regression (see ?sratio::scale_for_betareg)
scale_method <- "sv"

settings = list(
  n_sims = n_sims,
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
  n_cores = n_cores
)

# Function to run selectivity simulation

run_selectivity_simulation <- function(settings) {
  
  list2env(settings, envir = environment())
  
  # Draw samples
  catch_samples <- vector(mode = "list", length = n_sims)
  
  cat(paste0(as.POSIXct(Sys.time()), " Drawing catch samples\n"))
  
  for(ii in 1:n_sims) {
    
    paired_samples <- data.frame()
    
    for(jj in 1:n_pairs) {
      
      draw <- simulate_paired_sample(
        size = size, 
        abundance = abundance, 
        availability = availability, 
        demographic_comp_pars = demographic_comp_pars, 
        gear_q1 = gear_q1, 
        gear_q2 = gear_q2, 
        effort1 = effort1, 
        effort2 = effort2, 
        selectivity_opts1 = selectivity_opts1, 
        selectivity_opts2 = selectivity_opts2,
        return_vars = FALSE
      )$samples
      
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
    
    cat(paste0(as.POSIXct(Sys.time()), " Running analysis on sample ", kk, " out of ", n_sims, "\n"))
    
    sel_sim <- catch_samples[[kk]]
    
    sel_sim$size_bin <- round_any(sel_sim$size, accuracy = size_bin_width)
    
    # boot_samples[[kk]]  <- sratio::two_stage_bootstrap(
    #   count1 = sel_sim$n[sel_sim$treatment == 1],
    #   count2 = sel_sim$n[sel_sim$treatment == 2],
    #   size1 = sel_sim$size_bin[sel_sim$treatment == 1],
    #   size2 = sel_sim$size_bin[sel_sim$treatment == 2],
    #   block1 = sel_sim$haul[sel_sim$treatment == 1],
    #   block2 = sel_sim$haul[sel_sim$treatment == 2],
    #   n_draws = n_boot_samples,
    #   seed = NULL,
    #   treatment_name1 = 1,
    #   treatment_name2 = 2
    #   )
    
    cat(paste0(as.POSIXct(Sys.time()), " Drawing bootstrap samples\n"))
    
    boot_samples <- sratio::two_stage_bootstrap(count1 = sel_sim$n[sel_sim$treatment == 1],
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
    # for(ll in 1:length(boot_samples[[kk]])) {
    
    cat(paste0(as.POSIXct(Sys.time()), " Adding effort to bootstrap samples\n"))
    
    for(ll in 1:length(boot_samples)) {
      # boot_samples[[kk]][[ll]] <- dplyr::inner_join(boot_samples[[kk]][[ll]],
      #                                               data.frame(treatment = c(1, 2),
      #                                                          effort = c(effort1, effort2)),
      #                                               by = "treatment")
      
      boot_samples[[ll]] <- dplyr::inner_join(boot_samples[[ll]],
                                              data.frame(treatment = c(1, 2),
                                                         effort = c(effort1, effort2)),
                                              by = "treatment")
      
    }
    
    # Selectivity ratio model selection between binomial and beta binomial based on matchup level cross validation
    samples_wide <- catch_samples[[kk]] |>
      tidyr::pivot_wider(names_from = "treatment", 
                         values_from = "count", 
                         values_fill = 0, 
                         names_prefix = "count") |>
      dplyr::mutate(effort1 = effort1,
                    effort2 = effort2)
    
    cat(paste0(as.POSIXct(Sys.time()), " Stating model selection\n"))
    
    pratio <- sratio_cv(
      count1 = samples_wide$count1,
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
    

    cat(paste0(as.POSIXct(Sys.time()), " Fitting to bootstrap samples\n"))
    
    boot_sratio <- sratio::sratio_fit_bootstrap(
      x = boot_samples,
      # x = boot_samples[[kk]],
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
    
    cat(paste0(as.POSIXct(Sys.time()), " Summarising results\n"))
    
    boot_fit_sratio[[kk]] <- boot_sratio
    
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
  
  output <- list(
    bootstrap_results = boot_results,
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
                    n_cores = n_cores
    )
  )
  
  end_time <- Sys.time()
  
  cat(paste0("Total run time: ", round(difftime(end_time, start_time, units = "hour"), 1), " hrs"))
  
  return(output)
  
}

# Function to plot simulation results

plot_simulation_results <- function(x) {
  
  results_df <- do.call(what = "rbind", x$bootstrap_results)
  
  
  results_summary <- dplyr::mutate(results_df, 
                                resid = sratio_q500-true_sratio,
                                within_50 = true_sratio > sratio_q250 & true_sratio < sratio_q750,
                                within_95 = true_sratio > sratio_q025 & true_sratio < sratio_q975,
                                relative_precision_95 = (sratio_q975-sratio_q025)/true_sratio,
                                relative_precision_50 = (sratio_q750-sratio_q250)/true_sratio)
  
  plot_sim_sratio <- ggplot() +
    geom_path(data = results_df, 
              mapping = aes(x = size, y = sratio_q500, group = draw, color = "sample", alpha = "sample", linetype = "sample")) +
    geom_path(data = dplyr::group_by(results_summary, size) |>
                dplyr::summarise(mean = mean(sratio_q500, na.rm = TRUE),
                                 true = mean(true_sratio, na.rm = TRUE)) |>
                tidyr::pivot_longer(cols = c("mean", "true")),
              mapping = aes(x = size, y = value, color = name, alpha = name, linetype = name),
              linewidth = rel(1.1)) +
    scale_color_manual(values = c("mean" = "black",
                                  "true" = "red",
                                  "sample" = "grey50")) +
    scale_alpha_manual(values = c("mean" = 1,
                                  "true" = 1,
                                  "sample" = 0.25), guide = "none") +
    scale_linetype_manual(values = c("mean" = 1,
                                     "true" = 2,
                                     "sample" = 1), guide = "none") +
    scale_x_continuous(name = "Size") +
    scale_y_continuous(name = expression(S[21])) +
    theme_bw() +
    theme(legend.title = element_blank())
  
  
  plot_abs_bias <-  ggplot() +
    geom_path(data = dplyr::group_by(results_summary, size) |>
                dplyr::summarise(mean_bias = mean(resid, na.rm = TRUE)),
              mapping = aes(x = size, y = mean_bias),
              color = "red") +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(name = "Size") +
    scale_y_continuous(name = "Mean abs. bias") +
    theme_bw()
  
  
  plot_rel_bias <- ggplot() +
    geom_path(data = dplyr::group_by(results_summary, size, true_sratio) |>
                dplyr::summarise(mean_bias = mean(resid, na.rm = TRUE)) |>
                dplyr::mutate(rel_bias = mean_bias/true_sratio*100),
              mapping = aes(x = size, y = rel_bias),
              color = "red") +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_continuous(name = "Size") +
    scale_y_continuous(name = "Mean rel. bias (%)") +
    theme_bw()
  
  plot_within_ci <- ggplot() +
    geom_path(data = results_summary |>
                dplyr::group_by(size) |>
                dplyr::summarise(n = n(),
                                 within_50 = sum(within_50),
                                 within_95 = sum(within_95)) |>
                dplyr::mutate(`50%` = within_50/n*100,
                              `95%` = within_95/n*100) |>
                tidyr::pivot_longer(cols = c("50%", "95%")),
              mapping = aes(x = size, y = value, color = name)) +
    scale_y_continuous(name = "True within CI (%)", limits = c(0, 100)) +
    scale_x_continuous(name = "Size") +
    scale_color_colorblind(name = "CI") +
    theme_bw() +
    theme(legend.position = "none")
  
  plot_rel_precision <- ggplot() +
    geom_path(data = results_summary |>
                dplyr::group_by(size) |>
                dplyr::summarise(`95%` = mean(relative_precision_95),
                                 `50%` = mean(relative_precision_50)) |>
                tidyr::pivot_longer(cols = c("95%", "50%")),
              mapping = aes(x = size, y = value, color = name)) +
    scale_y_log10(name = "Rel. precision") +
    scale_x_continuous(name = "Size") +
    scale_color_colorblind(name = "CI") +
    theme_bw() +
    theme(legend.position = c(0.5, 0.8), legend.direction = "horizontal")
  
  # three_panel_bias_ci_prec <- cowplot::plot_grid(
  #   plot_rel_bias,
  #   plot_within_ci,
  #   plot_rel_precision,
  #   nrow = 3, align = "v"
  # )
  
  output <- list(plot_sim_sratio = plot_sim_sratio,
                 plot_abs_bias = plot_abs_bias,
                 plot_rel_bias = plot_rel_bias,
                 plot_within_ci = plot_within_ci,
                 plot_rel_precision = plot_rel_precision,
                 settings = x$settings
                 )
  
  return(output)
  
}

test_plots <- plot_simulation_results(x = test)





ggplot() +
  geom_path(data = results_df,
            mapping = aes(x = size, y = sratio_q500, group = draw), 
            color = "grey50", 
            alpha = 0.3) +
  geom_path(data = results_df |>
              dplyr::group_by(size) |>
              dplyr::summarise(mean_fit = mean(sratio_q500),
                               true = mean(true_sratio)) |>
              tidyr::pivot_longer(cols = c("mean_fit", "true")),
            mapping = aes(x = size, y = value, color = name, linetype = name), 
            linewidth = rel(1.1)) +
  scale_color_manual(values = c("true" = "red", "mean_fit" = "black")) +
  scale_y_continuous(name = expression(S[21])) +
  scale_x_continuous(name = "Size") +
  theme_bw() +
  theme(legend.title = element_blank())



