# Attempt to replicate Somerton et al. (2002)

# Somerton, D.A., Otto, R.S., Syrjala, S.E., 2002. Can changes in tow duration on bottom trawl surveys lead to changes in CPUE and mean size? Fish. Res. 55, 63–70. https://doi.org/10.1016/S0165-7836(01)00293-4

library(sratio)

crab_species_codes <- c(68560, 68580, 69322)

channel <- sratio::get_connected()

somerton_hauls_1998 <- 
  RODBC::sqlQuery(
    channel = channel,
    query = 
      paste0(
        "SELECT hauljoin, vessel, cruise, haul, duration, distance_fished, net_width 
      FROM racebase.haul
      WHERE
        cruise = 199801
        AND hauljoin IN (", paste(unique(sratio::crab_size_1995_1998$HAULJOIN), collapse = ", "),
        ")"
      ) 
  ) |>
  dplyr::mutate(
    AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
    TREATMENT = 
      factor(
        plyr::round_any(DURATION, 0.25), 
        levels = c(0.25, 0.5), 
        labels = c(15, 30)
      )
  ) |>
  dplyr::inner_join(
    dplyr::select(
      sratio::otto_key_1998, # Use tow pairs from Somerton et al. (2002)
      VESSEL, CRUISE, HAUL, TOW_PAIR, TOW_BLOCK
    ),
    by = c("VESSEL", "CRUISE", "HAUL")
  ) |>
  dplyr::arrange(TOW_BLOCK)

somerton_crab <- 
  somerton_hauls_1998 |>
  dplyr::select(HAULJOIN, VESSEL, CRUISE, HAUL, AREA_SWEPT_KM2, DURATION, TREATMENT, TOW_PAIR, TOW_BLOCK) |>
  dplyr::inner_join(
    dplyr::filter(sratio::crab_size_1995_1998),
    by = c("HAULJOIN", "VESSEL", "CRUISE", "HAUL")
    ) |>
      dplyr::filter(SPECIES_CODE %in% crab_species_codes,  # RKC, snow, Tanner from 1998
                    CRUISE == 199801) |>
      dplyr::mutate(SPECIES_CODE_SEX  = as.numeric(paste0(SPECIES_CODE, SEX)),  # Combine sex and species code
                    VESSEL = factor(VESSEL),
                    SEX = factor(SEX)) # Make vessel a factor for models

# Calculate weighted mean carapace width (snow and Tanner) or length (RKC)
mean_carapace_size <- 
  somerton_crab |>
  dplyr::group_by(
    HAULJOIN, 
    VESSEL, 
    CRUISE, 
    TREATMENT, 
    TOW_BLOCK, 
    HAUL, 
    SPECIES_CODE,
    SEX) |>
  dplyr::summarise(
    MEAN_WIDTH = sratio::weighted_mean(x = WIDTH, w = SAMPLING_FACTOR, na_rm = TRUE),
    MEAN_LENGTH = sratio::weighted_mean(x = LENGTH, w = SAMPLING_FACTOR, na_rm = TRUE),
    .groups = "keep"
    ) |>
  dplyr::mutate(MEAN_CARAPACE = ifelse(SPECIES_CODE == 69322, MEAN_LENGTH, MEAN_WIDTH))


# Carapace size models ----
mod_carapace_rkc_m <- 
  lm(
    formula = MEAN_CARAPACE ~ TREATMENT + VESSEL,
    data = mean_carapace_size[mean_carapace_size$SPECIES_CODE == 69322 & mean_carapace_size$SEX == 1, ]
  )

summary(mod_carapace_rkc_m)
anova(mod_carapace_rkc_m)

mod_carapace_rkc_f <- 
  lm(
  formula = MEAN_CARAPACE ~ TREATMENT + VESSEL,
  data = mean_carapace_size[mean_carapace_size$SPECIES_CODE == 69322 & mean_carapace_size$SEX == 2, ]
)

summary(mod_carapace_rkc_f)
anova(mod_carapace_rkc_f)

mod_carapace_tc_m <- 
  lm(
    formula = MEAN_CARAPACE ~ TREATMENT + VESSEL,
    data = mean_carapace_size[mean_carapace_size$SPECIES_CODE == 68560 & mean_carapace_size$SEX == 1, ]
  )

summary(mod_carapace_tc_m)
anova(mod_carapace_tc_m)

mod_carapace_tc_f <- 
  lm(
    formula = MEAN_CARAPACE ~ TREATMENT + VESSEL,
    data = mean_carapace_size[mean_carapace_size$SPECIES_CODE == 68560 & mean_carapace_size$SEX == 2, ]
  )

summary(mod_carapace_tc_f)
anova(mod_carapace_tc_f)

mod_carapace_sc_m <- 
  lm(
    formula = MEAN_CARAPACE ~ TREATMENT + VESSEL,
    data = mean_carapace_size[mean_carapace_size$SPECIES_CODE == 68580 & mean_carapace_size$SEX == 1, ]
  )

summary(mod_carapace_tc_m)
anova(mod_carapace_tc_m)

mod_carapace_sc_f <- 
  lm(
    formula = MEAN_CARAPACE ~ TREATMENT + VESSEL,
    data = mean_carapace_size[mean_carapace_size$SPECIES_CODE == 68580 & mean_carapace_size$SEX == 2, ]
  )

summary(mod_carapace_tc_f)
anova(mod_carapace_tc_f)

# CPUE ratio regression ----- 

cpue <- 
  somerton_crab |>
  dplyr::group_by(VESSEL, CRUISE, HAUL, TOW_PAIR, TREATMENT, AREA_SWEPT_KM2, SPECIES_CODE, SEX) |>
  dplyr::summarize(COUNT = sum(SAMPLING_FACTOR)) |>
  dplyr::ungroup() |>
  dplyr::mutate(CPUE_NO_KM2 = COUNT/AREA_SWEPT_KM2) |>
  dplyr::select(-HAUL, -CRUISE) |>
  tidyr::pivot_wider(names_from = "TREATMENT", 
                     values_from = c("CPUE_NO_KM2", "COUNT", "AREA_SWEPT_KM2"),
                     values_fill = 0) |>
  dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0) |># Remove pairs with zero CPUE
  dplyr::mutate(CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
                COMBINED_COUNT = COUNT_30 + COUNT_15)


# RKC models ----
mod_cpue_rkc_1 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX + VESSEL,
    data = cpue[cpue$SPECIES_CODE == 69322, ]
  )

summary(mod_cpue_rkc_1)
anova(mod_cpue_rkc_1)

mod_cpue_rkc_2 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX,
    data = cpue[cpue$SPECIES_CODE == 69322, ]
  )

mod_cpue_rkc_3 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ VESSEL,
    data = cpue[cpue$SPECIES_CODE == 69322, ]
  )

mod_cpue_rkc_4 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ 1,
    data = cpue[cpue$SPECIES_CODE == 69322, ]
  )

# Best model is #4
AIC(mod_cpue_rkc_1, mod_cpue_rkc_2, mod_cpue_rkc_3, mod_cpue_rkc_4)
summary(mod_cpue_rkc_4)

# Snow crab models ----
mod_cpue_sc_1 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX + VESSEL,
    data = cpue[cpue$SPECIES_CODE == 68580, ]
  )

summary(mod_cpue_sc_1)
anova(mod_cpue_sc_1)

mod_cpue_sc_2 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX,
    data = cpue[cpue$SPECIES_CODE == 68580, ]
  )

mod_cpue_sc_3 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ VESSEL,
    data = cpue[cpue$SPECIES_CODE == 68580, ]
  )

mod_cpue_sc_4 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ 1,
    data = cpue[cpue$SPECIES_CODE == 68580, ]
  )

# Best model is #4
AIC(mod_cpue_sc_1, mod_cpue_sc_2, mod_cpue_sc_3, mod_cpue_sc_4)
summary(mod_cpue_sc_4)


# Tanner crab models ----
mod_cpue_tc_1 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX + VESSEL,
    data = cpue[cpue$SPECIES_CODE == 68560, ]
  )

summary(mod_cpue_tc_1)
anova(mod_cpue_tc_1)

mod_cpue_tc_2 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX,
    data = cpue[cpue$SPECIES_CODE == 68560, ]
  )

mod_cpue_tc_3 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ VESSEL,
    data = cpue[cpue$SPECIES_CODE == 68560, ]
  )

mod_cpue_tc_4 <- 
  lm(
    formula = CPUE_LOG_RATIO ~ 1,
    data = cpue[cpue$SPECIES_CODE == 68560, ]
  )

# Best model is #4
AIC(mod_cpue_tc_1, mod_cpue_tc_2, mod_cpue_tc_3, mod_cpue_tc_4)
summary(mod_cpue_tc_4)

# Function to extract model intercept, variance, and bias-corrected ratio from Somerton's log-ratio models
extract_bias_corrected_ratio <- function(mod) {
  
  ratio_variance <- summary(mod)$sigma^2
  ratio_coef <- coef(mod)["(Intercept)"]
  ratio_estimate <- exp(ratio_coef + 0.5 * ratio_variance)
  
  return(
    c("ratio_intercept" = unname(ratio_coef),
      "ratio_variance" = unname(ratio_variance),
      "ratio_estimate" = unname(ratio_estimate))
  )
  
}

# The best model for all three species is an intercept-only model. Therefore, we use a log ratio 
# estimator to estimate the mean ratio.

extract_bias_corrected_ratio(mod_cpue_rkc_4)
extract_bias_corrected_ratio(mod_cpue_sc_4)
extract_bias_corrected_ratio(mod_cpue_tc_4)

log_ratio_estimator <- function(x, y, confint = 0.95) {
  
  stopifnot("x and y must be the same length." = length(x) == length(y))
  stopifnot("x and y must contain only positive values." = !any(x <= 0 | y <= 0))
  
  n <- length(x)
  
  log_mu <- mean(log(x/y), na.rm = TRUE)
  
  var_log_mu <- var(log(x/y), na.rm = TRUE)
  
  crit <- qt(1-(1-confint)/2, df = n-1)
  
  # Bias correction factor
  bc <- 0.5 * var_log_mu / n
  
  mu <- exp(log_mu + bc)
  
  se <- sqrt(var_log_mu)/sqrt(n)
  
  lci_mu <- exp(log_mu - se * crit + bc)
  
  uci_mu <- exp(log_mu + se * crit + bc)
  
  return(
    list(
      mu = mu, 
      lci_mu = lci_mu,
      uci_mu = uci_mu,
      log_mu = log_mu, 
      n = n, 
      se = se
    )
  )
  
}


log_ratio_estimator_bootstrap <- function(x, y, conf.level = 0.95, n_boot = 1000, return_boot = FALSE) {
  
  n <- length(x)
  log_mu <- log(x) - log(y)
  mean_log <- mean(log_mu)
  var_log <- var(log_mu)
  
  # Bias-corrected geometric mean ratio
  mu <- exp(mean_log + 0.5 * var_log / n)
  
  mean_var <- function(x) {
    cbind(mean(x), var(x))
  }
  
  # Bootstrap sampling
  boot_means <- replicate(n_boot, {
    indices <- sample(1:n, size = n, replace = TRUE)
    mean_var(log(x[indices]) - log(y[indices]))
  }, simplify = TRUE)
  
  boot_means <- t(boot_means)
  
  # Calculate bias correction factor for each sample
  boot_means <- cbind(boot_means, 0.5 * boot_means[, 2] / n)
  
  # Bootstrap confidence interval on the ratio scale
  ci_bounds <- quantile(
    exp(boot_means[, 1] + boot_means[, 3]), 
    probs = c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  )
  
  result <- list(
    mu = mu,
    lci_mu = ci_bounds[1],
    uci_mu = ci_bounds[2],
    log_mu = mean_log,
    n = n
  )
  
  if(return_boot) {
    result$bootstrap_distribution <- exp(boot_means)
  }
  
  return(result)
}


log_ratio_estimator_bootstrap(x = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68560],
                    y = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560],
              n_boot = 1000)

log_ratio_estimator(x = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68560],
                    y = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560])

log_ratio_estimator(x = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68580],
                    y = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580])

log_ratio_estimator(x = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 69322],
                    y = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322])


fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68580 & cpue$SEX == 1],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580 & cpue$SEX == 1],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68580 & cpue$SEX == 2],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580 & cpue$SEX == 2],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 69322 & cpue$SEX == 1],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322 & cpue$SEX == 1],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 69322 & cpue$SEX == 2],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322 & cpue$SEX == 2],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68560 & cpue$SEX == 1],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560 & cpue$SEX == 1],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68560 & cpue$SEX == 2],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560 & cpue$SEX == 2],
                 method = 4)





ggplot() +
  geom_point(data = cpue_wrong,
             mapping = aes(x = CPUE_NO_KM2_30, y = CPUE_NO_KM2_15)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~SPECIES_CODE)


cowplot::plot_grid(
  ggplot() +
    geom_histogram(data = cpue,
                   mapping = aes(x = CPUE_NO_KM2_15/CPUE_NO_KM2_30)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~SPECIES_CODE),
  ggplot() +
    geom_histogram(data = cpue_wrong,
                   mapping = aes(x = CPUE_NO_KM2_15/CPUE_NO_KM2_30)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~SPECIES_CODE),
  nrow = 2
)


ggplot() +
  stat_ecdf(data = cpue,
                 mapping = aes(x = CPUE_NO_KM2_15/CPUE_NO_KM2_30)) +
  # scale_x_log10() +
  # scale_y_log10() +
  facet_wrap(~SPECIES_CODE)
