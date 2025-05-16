# Fishing Power Comparisons on numerical CPUE

library(ggplot2)
library(sratio)
library(cowplot)
library(fishmethods) # Version 1.13-1

spp_code = 685801

# Format data ----
dat <- 
  readRDS(file = here::here("analysis", "15_30", "output", "catch_at_length_1530.rds")) |>
  dplyr::mutate(TOTAL_FREQUENCY = SAMPLING_FACTOR * FREQUENCY) |>
  dplyr::group_by(HAULJOIN, SPECIES_CODE, MATCHUP, TREATMENT, AREA_SWEPT_KM2) |>
  dplyr::summarise(TOTAL_COUNT = sum(TOTAL_FREQUENCY), .groups = "keep") |>
  dplyr::ungroup() |>
  dplyr::mutate(CPUE_NO_KM2 = TOTAL_COUNT/AREA_SWEPT_KM2) |>
  dplyr::select(-HAULJOIN) |>
  tidyr::pivot_wider(
    values_from = c("AREA_SWEPT_KM2", "TOTAL_COUNT", "CPUE_NO_KM2"),
    names_from = "TREATMENT"
  ) |>
  dplyr::mutate(
    LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
    LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
    CPUE_RATIO = CPUE_NO_KM2_30/CPUE_NO_KM2_15,
    CPUE_LOG_RATIO = log(CPUE_RATIO)
  )


sel_dat <- dat |> 
  dplyr::filter(
  SPECIES_CODE == spp_code
)

# Standard log ratio estimator ----

set.seed(1337)

num_fpc_log_boot <- 
  log_ratio_estimator_bootstrap(
  x = sel_dat$CPUE_NO_KM2_30,
  y = sel_dat$CPUE_NO_KM2_15,
  n_boot = 1000,
  conf.level = 0.95,
  return_boot = TRUE
)


# Somerton ratio estimator ---- 
# No vessel or sex effects; Miller et al. (1984) bias correction

# Function to extract model intercept, variance, and bias-corrected ratio from Somerton's log-ratio models
somerton_bias_correction <- function(mod) {
  
  ratio_variance <- summary(mod)$sigma^2
  ratio_coef <- coef(mod)["(Intercept)"]
  ratio_estimate <- exp(ratio_coef + 0.5 * ratio_variance)
  
  return(
    c("ratio_intercept" = unname(ratio_coef),
      "ratio_variance" = unname(ratio_variance),
      "ratio_estimate" = unname(ratio_estimate))
  )
  
}

mod_somerton <- 
  lm(
    formula = CPUE_LOG_RATIO ~ 1,
    data = sel_dat
  )

num_fpc_somerton <- somerton_bias_correction(mod_somerton)


# Kappenman scale estimator and randomized block ANOVA ----

num_fpc_kappenman <-
  fishmethods::fpc(
    cpue1 = sel_dat$CPUE_NO_KM2_15,
    cpue2 = sel_dat$CPUE_NO_KM2_30,
    boot_type = "paired",
    decimals = 6,
    nboot = 1000,
    method = c(2, 4)
  )


# Bayesian zero-intercept linear regression ----

num_mod_zeroint <- 
    brms::brm(
      formula = LOG_CPUE_NO_KM2_30 ~ LOG_CPUE_NO_KM2_15 + 0, 
      data = sel_dat, 
      iter = 10000, 
      chains = 4,
      thin = 5,
      warmup = 2000
    )
  
  posterior_df <- brms::as_draws_df(num_mod_zeroint)
  
  posterior_df$SPECIES_CODE <- spp_code
  
  # Diagnostics
  pp_check(num_mod_zeroint, type = "dens_overlay", ndraws = 1000)
  pp_check(num_mod_zeroint, type = "loo_pit_overlay", ndraws = 100)
  pp_check(num_mod_zeroint, type = "dens_overlay", ndraws = 1000)
  pp_check(num_mod_zeroint, type = "scatter_avg", ndraws = 100)
  pp_check(num_mod_zeroint, type = "loo_pit_qq", ndraws = 4000, moment_match = TRUE)
  pp_check(num_mod_zeroint, type = "pit_ecdf", ndraws = 4000)
  
  loo_check <- brms::loo(num_mod_zeroint, moment_match = TRUE)
  
  loo::pareto_k_values(loo_check)
  
  # Generate vector of log(CPUE[15])
  new_x <- data.frame(
    LOG_CPUE_NO_KM2_15 = 
      seq(
        plyr::round_any(
          min(sel_dat$LOG_CPUE_NO_KM2_15), 
          accuracy = 0.1, 
          f = floor
        ), 
        plyr::round_any(
          max(sel_dat$LOG_CPUE_NO_KM2_15), 
          accuracy = 0.1, 
          f = ceiling), 
        by = 0.01)
  )
  
  new_x$CPUE_NO_KM2_15 <- exp(new_x$LOG_CPUE_NO_KM2_15)
  
  # Get posterior linear predictions (log scale, without random effects)
  fit_vals <- 
    posterior_linpred(
      num_mod_zeroint, 
      newdata = new_x, 
      transform = FALSE, 
      re_formula = NA
    )
  
  # Compute mu (mean of log preds) and sigma^2 (variance) for each x
  mu <- colMeans(fit_vals)
  sigma2 <- apply(fit_vals, 2, var)
  
  # Bias-corrected mean predictions
  bias_corrected <- exp(mu + 0.5 * sigma2)
  
  # Bias-corrected credible intervals for predictions
  lower <- apply(exp(fit_vals + 0.5 * sigma2), 2, quantile, probs = 0.025)
  upper <- apply(exp(fit_vals + 0.5 * sigma2), 2, quantile, probs = 0.975)
  
  bc_preds <- new_x |>
    dplyr::mutate(
      SPECIES_CODE = spp_code,
      CPUE_NO_KM2_30_FIT = bias_corrected,
      CPUE_NO_KM2_30_LWR = lower,
      CPUE_NO_KM2_30_UPR = upper
    )
  
  p1 <- ggplot() +
    geom_ribbon(data = bc_preds,
                mapping = aes(x = CPUE_NO_KM2_15, ymin = CPUE_NO_KM2_30_LWR, ymax = CPUE_NO_KM2_30_UPR), alpha = 0.5) +
    geom_path(data = bc_preds,
              mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_30_FIT)) +
    geom_point(data = sel_dat,
               mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_30)) +
    geom_abline(slope = 1, intercept = 0, linetype = 3) +
    geom_abline(slope = c(1.2, 0.8), intercept = 0, linetype = 3) +
    ggtitle(label = sratio::species_code_label(spp_code, type = "common_name")) +
    scale_x_log10(name = expression(CPUE[15]*' '*('#'%.%km^-2))) +
    scale_y_log10(name = expression(CPUE[30]*' '*('#'%.%km^-2))) +
    theme_bw()
  
  print(p1)

# VAST ----