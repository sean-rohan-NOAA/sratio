library(sratio)
library(fitdistrplus)

# Successful slope tows

slope_tows_83112 <- dplyr::filter(sratio::data_ss$haul, GEAR == 44, BOTTOM_DEPTH > 180)
nrow(slope_tows_83112)

sratio::data_ss$haul |> dplyr::filter(GEAR == 172)

# Gear codes and species names ----
gear_codes <- 
  data.frame(
    GEAR = c("PNE", "83-112"),
    TREATMENT = factor(c(172, 44))
  )

species_codes <-
  data.frame(
    SPECIES_CODE = c(21740, 21720, 20510, 30060, 30420, 30152, 30051, 30052, 30020, 10115, 10110, 10112, 10130, 471, 420, 435, 440, 455, 472, 475, 477, 480, 485,  68580),
    COMMON_NAME = c("walleye pollock", "Pacific cod", "sablefish", "Pacific ocean perch", "northern rockfish", "dusky rockfish", "BS/RE rockfish",  "BS/RE rockfish", "shortspine thornyhead", "Greenland turbot", "arrowtooth flounder", "Kamchatka flounder", "flathead sole", rep("skates", 10),"snow crab")
  )

species_codes$COMMON_NAME <- 
  factor(
    species_codes$COMMON_NAME, 
    levels = species_codes$COMMON_NAME, 
    labels = species_codes$COMMON_NAME
  )

# set analysis species based on sample sizes
analysis_species <- c(
  "walleye pollock", 
  "Pacific cod", 
  "sablefish", 
  "Pacific ocean perch", 
  "BS/RE rockfish",  
  "shortspine thornyhead", 
  "Greenland turbot", 
  "arrowtooth flounder", 
  "Kamchatka flounder", 
  "flathead sole", 
  "skates", 
  "snow crab"
)

# Get start and end dates for hauls in each year
(fishing_dates <- sratio::data_ss$haul |> 
  dplyr::group_by(YEAR, GEAR, HAUL_TYPE) |>
  dplyr::summarise(min_date = min(START_TIME),
                   max_date = max(START_TIME),
                   n = n()))

# Format CPUE data for model fitting ----

cpue_data <- 
  data_ss$catch |>
  dplyr::inner_join(
    species_codes
  ) |>
  dplyr::select(-SPECIES_CODE) |>
  dplyr::group_by(VESSEL, CRUISE, HAUL, HAULJOIN, MATCHUP, COMMON_NAME) |>
  dplyr::summarise(WEIGHT = sum(WEIGHT),
                   NUMBER_FISH = sum(NUMBER_FISH)) |>
  dplyr::ungroup() |>
  dplyr::inner_join(
    dplyr::select(
      data_ss$haul, VESSEL, CRUISE, HAUL, MATCHUP, TREATMENT, AREA_SWEPT_KM2
    )
  ) |>
  dplyr::inner_join(gear_codes) |>
  dplyr::mutate(
    CPUE_KGKM2 = WEIGHT/AREA_SWEPT_KM2,
    CPUE_NOKM2 = NUMBER_FISH/AREA_SWEPT_KM2,
) |>
  dplyr::select(-HAULJOIN, -VESSEL, -HAUL, -CRUISE)

cpue_wide <- cpue_data |>
  dplyr::select(-GEAR) |>
  tidyr::pivot_wider(
    values_from = c("CPUE_KGKM2", "CPUE_NOKM2", "NUMBER_FISH", "WEIGHT", "AREA_SWEPT_KM2"),  
    names_from = TREATMENT, 
    values_fill = 0
  ) |>
  dplyr::filter(COMMON_NAME %in% analysis_species) |>
  dplyr::mutate(CPUE_LOG_RATIO = log(CPUE_KGKM2_172/CPUE_KGKM2_44),
                log_44 = log(CPUE_KGKM2_44+1),
                log_172 = log(CPUE_KGKM2_172+1))


write.csv(cpue_wide, file = here::here("analysis", "shelf_slope", "data", "shelf_slope_cpue.csv"), row.names = FALSE)

# Draw bootstrap samples for each species
bootstrap_matchups <- vector(mode = "list", length = 1000)
set.seed(999)
for(ii in 1:1000) {
  
  bootstrap_matchups[[ii]] <- 
    dplyr::select(cpue_wide, COMMON_NAME, MATCHUP) |>
    dplyr::group_by(COMMON_NAME) |>
    dplyr::reframe(n = n(), 
                   MATCHUP = sample(MATCHUP, size = n, replace = TRUE)) |>
    dplyr::mutate(draw = ii) |>
    dplyr::select(-n) |>
    dplyr::left_join(cpue_wide, join_by(COMMON_NAME, MATCHUP))
  
}

saveRDS(bootstrap_matchups, file = here::here("analysis", "shelf_slope", "output", "fpc_bootstrap_samples.rds"))



# Haul summary table ----
haul_geometry <- 
  sratio::data_ss$haul |> 
  dplyr::group_by(TREATMENT) |>
  dplyr::summarise(MEAN_SPREAD = format(round(mean(NET_WIDTH), 1), nsmall = 1),
                   MEAN_HEIGHT = format(round(mean(NET_HEIGHT), 1), nsmall = 1),
                   MIN_SPREAD = format(round(min(NET_WIDTH), 1), nsmall = 1),
                   MAX_SPREAD = format(round(max(NET_WIDTH), 1), nsmall = 1),
                   MIN_HEIGHT = format(round(min(NET_HEIGHT), 1), nsmall = 1),
                   MAX_HEIGHT = format(round(max(NET_HEIGHT), 1), nsmall = 1),
                   MEAN_AREA_SWEPT_KM2 = format(round(mean(AREA_SWEPT_KM2), 3), nsmall = 3),
                   MIN_AREA_SWEPT_KM2 = format(round(min(AREA_SWEPT_KM2), 3), nsmall = 3),
                   MAX_AREA_SWEPT_KM2 = format(round(max(AREA_SWEPT_KM2), 3), nsmall = 3),
                   MEAN_DURATION = format(round(mean(DURATION*60), 1), nsmall = 1),
                   MIN_DURATION = format(round(min(DURATION*60), 1), nsmall = 1),
                   MAX_DURATION = format(round(max(DURATION*60), 1), nsmall = 1),
                   MEAN_SPEED = format(round(mean(DISTANCE_FISHED/DURATION)/1.852, 1), nsmall = 1),
                   MIN_SPEED = format(round(min(DISTANCE_FISHED/DURATION)/1.852, 1), nsmall = 1),
                   MAX_SPEED = format(round(max(DISTANCE_FISHED/DURATION)/1.852, 1), nsmall = 1),
                   n = n()) |>
  dplyr::inner_join(gear_codes) |>
  dplyr::mutate(
    Spread = paste0(MEAN_SPREAD, " (", MIN_SPREAD, "-", MAX_SPREAD, ")"),
    Height = paste0(MEAN_HEIGHT, " (", MIN_HEIGHT, "-", MAX_HEIGHT, ")"),
    `Area swept` = paste0(MEAN_AREA_SWEPT_KM2, " (", MIN_AREA_SWEPT_KM2, "-", MAX_AREA_SWEPT_KM2, ")"),
    Duration = paste0(MEAN_DURATION, " (", MIN_DURATION, "-", MAX_DURATION, ")"),
    `Tow speed` = paste0(MEAN_SPEED, " (", MIN_SPEED, "-", MAX_SPEED, ")")
  ) |>
  dplyr::select(
    Gear = GEAR,
    Spread,
    Height,
    `Area swept`,
    Duration,
    `Tow speed`
  )

write.csv(
  haul_geometry, 
  file = here::here("analysis", "shelf_slope", "plots", "haul_geometry.csv"),
  row.names = FALSE
)


# Size summary table ----

size_summary <-
  data_ss$size |>
  dplyr::inner_join(species_codes) |>
  dplyr::select(-SPECIES_CODE) |>
  dplyr::inner_join(dplyr::select(data_ss$haul, MATCHUP, HAULJOIN, TREATMENT)) |>
  dplyr::inner_join(gear_codes) |>
  dplyr::mutate(SIZE = ifelse(is.na(LENGTH), WIDTH, LENGTH)) |>
  dplyr::group_by(
    COMMON_NAME, GEAR
  ) |>
  dplyr::summarise(
    n_lengths = n(),
    MEAN_SIZE = format(round(mean(SIZE, na.rm = TRUE), 1), nsmall = 1),
    MIN_SIZE = format(round(min(SIZE, na.rm = TRUE), 1), nsmall = 1),
    MAX_SIZE = format(round(max(SIZE, na.rm = TRUE), 1), nsmall = 1)
  ) |>
  dplyr::mutate(SIZE = paste0(MEAN_SIZE, " (",MIN_SIZE, "-", MAX_SIZE, ")")) |>
  dplyr::select(-MEAN_SIZE, -MIN_SIZE, -MAX_SIZE)

species_hauls <- 
  cpue_data |>
  dplyr::group_by(COMMON_NAME, GEAR) |>
  dplyr::summarise(
    n_positive = sum(CPUE_NOKM2 > 0)
  ) |> 
  dplyr::mutate(
    p_positive = n_positive/40 * 100
  ) |>
  dplyr::select(-p_positive)

cpue_hauls <- 
  cpue_data |>
  dplyr::group_by(COMMON_NAME, GEAR) |>
  dplyr::summarise(
    MEAN_CPUE_MTKM2 = format(round(mean(CPUE_KGKM2/1000, na.rm = TRUE), 3), nsmall = 3),
    MIN_CPUE_MTKM2 = format(round(min(CPUE_KGKM2/1000, na.rm = TRUE), 3), nsmall = 3),
    MAX_CPUE_MTKM2 = format(round(max(CPUE_KGKM2/1000, na.rm = TRUE), 3), nsmall = 3),
    MEAN_CPUE_NOKM2 = format(round(mean(CPUE_NOKM2/1000, na.rm = TRUE), 3), nsmall = 3),
    MIN_CPUE_NOKM2 = format(round(min(CPUE_NOKM2/1000, na.rm = TRUE), 3), nsmall = 3),
    MAX_CPUE_NOKM2 = format(round(max(CPUE_NOKM2/1000, na.rm = TRUE), 3), nsmall = 3)
  ) |>
  dplyr::mutate(
    CPUE_MTKM2 = paste0(MEAN_CPUE_MTKM2, "\n(", MIN_CPUE_MTKM2, "-", MAX_CPUE_MTKM2, ")"),
    CPUE_NOKM2 = paste0(MEAN_CPUE_NOKM2, "\n(", MIN_CPUE_NOKM2, "-", MAX_CPUE_NOKM2, ")")
    ) |>
  dplyr::select(-MEAN_CPUE_MTKM2, -MIN_CPUE_MTKM2, -MAX_CPUE_MTKM2, -MEAN_CPUE_NOKM2, -MIN_CPUE_NOKM2, -MAX_CPUE_NOKM2)
  

summary_table <- 
  dplyr::inner_join(species_hauls, cpue_hauls) |>
  dplyr::inner_join(size_summary)

write.csv(
  summary_table, 
  file = here::here("analysis", "shelf_slope", "plots", "catch_summary_table.csv"), 
  row.names = FALSE
)


# Map of samples ----

akgfmaps::get_base_layers(select.region = c("sebs", "ebs.slope"), set.crs = "EPSG:3338")



# Setup functions for each method ----

# Fit GAMs models ----


fit_mgcv <- function(x, common_name, formula = log_172 ~ s(log_44, bs = "tp"), family = gaussian(), loocv = FALSE) {
  
  x_sel <- dplyr::filter(x, COMMON_NAME == common_name)
  
  mod <- mgcv::gam(formula = formula, family = family, data = x_sel)
  
  predictor_name <- all.vars(formula)[2]

  newdata <- 
    data.frame(
      COMMON_NAME = common_name,
      pred = seq(min(x_sel[[predictor_name]]), max(x_sel[[predictor_name]]), 0.05)
    )
  
  names(newdata)[2] <- predictor_name
  
  fit <- predict(mod, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$fit <- fit$fit
  newdata$se.fit <- fit$se.fit
  newdata$lwr <- fit$fit - 2*fit$se.fit
  newdata$upr <- fit$fit + 2*fit$se.fit
  
  # Leave one out cross validation
  loocv_results <- NULL
  
  if(loocv) {
    
    loocv_results <- data.frame()
    
    for(ii in 1:nrow(x_sel)) {
      
      validation <- x_sel[ii, ]
      fitting <- x_sel[-ii, ]
      
      mod_loocv <- mgcv::gam(formula = formula, family = family, data = fitting)
      fold_fit <- predict(mod_loocv, newdata = validation, type = "response", se.fit = TRUE)
      validation$fit <- fold_fit$fit
      validation$se.fit <- fold_fit$se.fit
      validation$lwr <- fold_fit$fit - 2*fold_fit$se.fit
      validation$upr <- fold_fit$fit + 2*fold_fit$se.fit
      
      loocv_results <- rbind(loocv_results, validation)
      
    }
    
  }
  
  return(
    list(model = mod,
         loocv_results = loocv_results,
         fit = newdata)
  )
  
}

# Fishing power comparison methods ----
# Methods: 1 = Ratio of Means, 2 = Randomized Block ANOVA, 3 = Multiplicative Model, 4 = Kappenman 1992.

est_fpc <- function(x, common_name, cpue1_name = "CPUE_KGKM2_172", cpue2_name = "CPUE_KGKM2_44", method = 1:4, loocv = FALSE, kapp_zeros = "ind", boot_type = "paired", nboot = 1000) {
  
  x_sel <- dplyr::filter(x, COMMON_NAME == common_name)
  
  fpc_results <- 
    fishmethods::fpc(
      cpue1 = x_sel[[cpue1_name]], 
      cpue2 = x_sel[[cpue2_name]], 
      method = method, 
      nboot = nboot
    )
  
  fpc_results$COMMON_NAME <- common_name
  
  # Leave one out cross validation
  loocv_results <- NULL
  
  if(loocv) {
    
    loocv_results <- data.frame()
    
    for(ii in 1:nrow(x_sel)) {
      
      print(Sys.time())
      
      fold_fit <-
        cbind(x_sel[ii, ], 
              fishmethods::fpc(
                cpue1 = x_sel[[cpue1_name]][-ii], 
                cpue2 = x_sel[[cpue2_name]][-ii], 
                method = method,
                boot_type = boot_type,
                nboot = nboot
                )
        )
      
      loocv_results <- rbind(loocv_results, fold_fit)
      
    }
    
  }
  
  return(list(
    fpc = fpc_results,
    loocv_results = loocv_results)
  )
  
}


# Fit Somerton models ----

fit_somerton_ratio <- function(x, common_name, cpue1_name, cpue2_name, loocv = FALSE) {
  
  x_sel <- dplyr::filter(x, COMMON_NAME == common_name)
  
  ratio_mod <- 
    lm(
      formula = CPUE_LOG_RATIO ~ 1,
      data = x_sel
    )
  
  ratio <- miller_bias_correct(ratio_mod)
  
  ratio$COMMON_NAME <- common_name
  
  # Leave one out cross validation
  loocv_results <- NULL
  
  if(loocv) {
    
    loocv_results <- data.frame()
    
    for(ii in 1:nrow(x_sel)) {
      
      fold_mod <- 
        lm(
          formula = CPUE_LOG_RATIO ~ 1,
          data = x_sel[-ii,]
        )
      
      fold_fit <- 
        cbind(
          x_sel[ii,],
          miller_bias_correct(fold_mod)
        )
      
      loocv_results <- rbind(loocv_results, fold_fit)
      
    }
    
  }
  
  return(
    list(model = ratio_mod,
         ratio = ratio,
         loocv_results = loocv_results)
  )
  
}

# Function to extract model intercept, variance, and bias-corrected ratio from Somerton log-ratio models
miller_bias_correct <- function(mod) {
  
  if(all(names(coef(mod)) == "(Intercept)")) {
    
    log_ratio <- coef(mod)["(Intercept)"]
    var <- summary(mod)$sigma^2
    se <- summary(mod)$coefficients[, 2]
    ratio <- exp(log_ratio)
    ratio_bc <- exp(log_ratio + 0.5 * var)
    ratio_lci <- exp(log_ratio - 2 * se + 0.5 * var)
    ratio_uci <- exp(log_ratio + 2 * se + 0.5 * var)
    
    output <- data.frame(
      log_ratio = unname(log_ratio),
      var = unname(var),
      ratio = unname(ratio),
      ratio_bc = unname(ratio_bc),
      ratio_lci = ratio_lci,
      ratio_uci = ratio_uci
    )
    
    rownames(output) <- NULL
    
  } else {
    
    log_ratio <- coef(mod)[[1]]
    
    fit <- data.frame(
      LOG_CPUE_NO_KM2_30 = 
        seq(min(mod$model[,2]), max(mod$model[,2]), by = 0.01)
    )
    
    fit$log_ratio <- predict(mod, newdata = fit)
    fit$se <- predict(mod, newdata = fit, se.fit = TRUE)$se.fit
    var <- summary(mod)$sigma^2
    fit$fit <- exp(fit$log_ratio)
    fit$fit_bc <- exp(fit$log_ratio + 0.5 * var)
    fit$ratio_lci <- exp(fit$log_ratio - 2 * fit$se + 0.5 * var)
    fit$ratio_uci <- exp(fit$log_ratio + 2 * fit$se + 0.5 * var)
    
    output <-
      list(fit = fit,
           pars = c("var" = var, "log_ratio" = log_ratio))
    
  }
  
  return(output)
  
}

# Catch comparison rate ---- 
# Sensu O'Leary et al. (2021); analogous to catch comparison rate 

est_ccr <- function(x, common_name, cpue1_name = "CPUE_KGKM2_172", cpue2_name = "CPUE_KGKM2_44", loocv = FALSE, bootstrap_ci = 1000) {
  
  x_sel <- dplyr::filter(x, COMMON_NAME == common_name)
  
  n_obs <- nrow(x_sel)
  
  f_ccr <- function(c1, c2, fn = mean) {
    return(fn(c1/(c1+c2), na.rm = TRUE))
  }
  
  ccr_to_fpc <- function(val) {
    return(val/(1-val))
  }
  
  ccr <- 
    data.frame(
      COMMON_NAME = common_name,
      ccr_mean = f_ccr(c1 = x_sel[[cpue1_name]], c2 = x_sel[[cpue2_name]], fn = mean),
      ccr_sd = f_ccr(c1 = x_sel[[cpue1_name]], c2 = x_sel[[cpue2_name]], fn = sd),
      lwr_ci_boot = NA,
      upr_ci_boot = NA,
      ccr_median = f_ccr(c1 = x_sel[[cpue1_name]], c2 = x_sel[[cpue2_name]], fn = median)
    )
  
  ccr$fpc <- ccr_to_fpc(ccr$ccr_mean)
  ccr$fpc_median <- ccr_to_fpc(ccr$ccr_median)

  # bootstrap for CIs when bootstrap_ci is numeric
  
  if(is.numeric(bootstrap_ci)) {
    
    ccr_samp <- numeric(length = bootstrap_ci)
    
    for(jj in 1:bootstrap_ci) {
      
      boot_dat <- x_sel[sample(1:n_obs, size = n_obs, replace = TRUE), c(cpue1_name, cpue2_name)]
      
      ccr_samp[jj] <- f_ccr(c1 = boot_dat[[cpue1_name]], c2 = boot_dat[[cpue2_name]])
      
    }
    
    ccr$ccr_lwr_ci_boot <- quantile(ccr_samp, p = 0.025)
    ccr$ccr_upr_ci_boot <- quantile(ccr_samp, p = 0.975)
    
    ccr$fpc_lwr_ci_boot <- ccr_to_fpc(ccr$ccr_lwr_ci_boot)
    ccr$fpc_upr_ci_boot <- ccr_to_fpc(ccr$ccr_upr_ci_boot)
    
  }
  
  # Leave one out cross validation
  loocv_results <- NULL
  
  if(loocv) {
    
    loocv_results <- data.frame()
    
    for(ii in 1:nrow(x_sel)) {
      
      ccr_fit <- f_ccr(c1 = x_sel[[cpue1_name]][-ii], c2 = x_sel[[cpue2_name]][-ii])
      
      loocv_results <- rbind(loocv_results, cbind(x_sel[ii, ], ccr_fit))
      
    }
    
    loocv_results$fpc_fit <- ccr_to_fpc(loocv_results$ccr_fit)
    
  }
  
  return(
    list(
      ccr = ccr,
      loocv_results = loocv_results
    )
  )
  
  
}


# Median ratio (O'Leary  et al., 2021) ----

est_median_ratio <- function(x, common_name, cpue1_name, cpue2_name, loocv = FALSE, bootstrap_ci = 1000) {
  
  x_sel <- dplyr::filter(x, COMMON_NAME == common_name)
  
  x_cv <- x_sel  
  
  # Only positive values
  
  x_sel <- x_sel[x_sel[[cpue1_name]] > 0 & x_sel[[cpue2_name]] > 0, ]
  
  n_obs <- nrow(x_sel)
  
  ser <- 
    data.frame(
      COMMON_NAME = common_name,
      ser_median = median(x_sel[[cpue1_name]]/x_sel[[cpue2_name]], na.rm = TRUE),
      ser_mean = mean(x_sel[[cpue1_name]]/x_sel[[cpue2_name]], na.rm = TRUE),
      ser_lwr_ci_boot = NA,
      ser_upr_ci_boot = NA,
      n_pair = nrow(x_sel)
    )
  
  loocv_results <- NULL
  
  # Bootstrap to get CIs
  if(is.numeric(bootstrap_ci)) {
    
    ser_samp <- numeric(length = bootstrap_ci)
    
    for(jj in 1:bootstrap_ci) {
      
      boot_dat <- x_sel[sample(1:n_obs, size = n_obs, replace = TRUE), c(cpue1_name, cpue2_name)]
      
      ser_samp[jj] <- median(boot_dat[[cpue1_name]]/boot_dat[[cpue2_name]], na.rm = TRUE)
      
    }
    
    ser$ser_lwr_ci_boot <- quantile(ser_samp, p = 0.025)
    ser$ser_upr_ci_boot <- quantile(ser_samp, p = 0.975)
    
  }
  
  
  # Leave one out cross validation
  loocv_results <- NULL
  
  if(loocv) {
    
    loocv_results <- data.frame()
    
    for(ii in 1:nrow(x_cv)) {
      
      x_cv_fold <- x_cv[-ii, ]
      
      x_cv_fold <- x_cv_fold[x_cv_fold[[cpue1_name]] > 0 & x_cv_fold[[cpue2_name]] > 0, ]
      
      ser_median <- median(x_cv_fold[[cpue1_name]]/x_cv_fold[[cpue2_name]], na.rm = TRUE)
      
      loocv_results <- rbind(loocv_results, cbind(x_cv[ii, ], ser_median))
      
    }
    
  }
  
  return(
    list(
      ser = ser,
      loocv_results = loocv_results
    )
  )
  
}

# Catch comparison rate using beta distribution ----

est_ccr_beta <- function(x, common_name, cpue1_name = "CPUE_KGKM2_172", cpue2_name = "CPUE_KGKM2_44", loocv = FALSE, bootstrap_ci = 1000) {
  
  x_sel <- dplyr::filter(x, COMMON_NAME == common_name)
  
  n_obs <- nrow(x_sel)
  
  f_ccr <- function(c1, c2) {
    return(c1/(c1+c2))
  }
  
  ccr_to_fpc <- function(val) {
    return(val/(1-val))
  }
  
  # Catch comparison rate
  p <- f_ccr(c1 = x_sel[[cpue1_name]], c2 = x_sel[[cpue2_name]])
  
  beta_fit <- fitdistrplus::fitdist(p, "beta")
  
  ccr <- 
    data.frame(
      COMMON_NAME = common_name,
      par_alpha = beta_fit$estimate[1],
      par_beta = beta_fit$estimate[2]
    )
  
  ccr$p <- ccr$par_alpha/(ccr$par_alpha+ccr$par_beta)
  
  ccr$fpc <- ccr_to_fpc(ccr$p)
  
  
  # bootstrap for CIs when bootstrap_ci is numeric
  
  if(is.numeric(bootstrap_ci)) {
    
    boot_fpc <- numeric(length = bootstrap_ci)
    
    for(jj in 1:bootstrap_ci) {
      
      boot_p <- sample(p, size = n_obs, replace = TRUE)
      
      boot_fit <- fitdistrplus::fitdist(boot_p, "beta")
      
      boot_alpha <- boot_fit$estimate[1]
      
      boot_beta <- boot_fit$estimate[2]
      
      boot_p <- boot_alpha/(boot_alpha+boot_beta)
      
      boot_fpc[jj] <- ccr_to_fpc(boot_p)
      
    }
    
    ccr$fpc_lwr_ci_boot <- quantile(boot_fpc, p = 0.025)
    ccr$fpc_upr_ci_boot <- quantile(boot_fpc, p = 0.975)
    
  }
  
  # Leave one out cross validation
  loocv_results <- NULL
  
  if(loocv) {
    
    loocv_results <- data.frame()
    
    for(ii in 1:nrow(x_sel)) {
      
      loocv_fit <- fitdistrplus::fitdist(p[-ii], "beta")
      loocv_alpha <- loocv_fit$estimate[1]
      loocv_beta <- loocv_fit$estimate[2]
      fpc_fit <- loocv_alpha/(loocv_alpha+loocv_beta)
      
      loocv_results <- rbind(loocv_results, cbind(x_sel[ii, ], fpc_fit))
      
    }
    
  }
  
  return(
    list(
      ccr = ccr,
      loocv_results = loocv_results
    )
  )
  
  
}

# Loop through species to estimate relative fishing powers for biomass  ----

# Methods used:
# GAM, median ratio, log ratio with bias correction (Somerton), RCI/CCR beta (Fanning), 
# RCI/CCR mean (O'Leary), ratio median (O'Leary), Kappenman

# For each method:
# - Calculate estimator
# - Bootstrap to estimate CIs
# - Generate LOOCV predictions for method comparison

bootstrap_species <- 
  do.call(rbind, bootstrap_matchups)

gam_fits <- fpc_fits <- somerton_fits <- fpc_fits <- ccr_fits <- ccr_beta_fits <- 
  median_ratio_fits <- vector(mode = "list", length = length(analysis_species))

bootstrap_fits <- 

gam_pred <- data.frame()

loocv_pred <- data.frame()

set.seed(999)
for(ii in 1:1) { # Takes ~2 hours to run because of the Kappenman estimator; otherwise less than a minute
# for(ii in 1:length(analysis_species)) { # Takes ~2 hours to run because of the Kappenman estimator; otherwise less than a minute
  
  cat(analysis_species[ii], "\n")
  
  boot_cpue_samples <- 
    lapply(bootstrap_matchups, dplyr::filter, COMMON_NAME == analysis_species[[ii]])
  
  # GAM - No specific author, but see (Jones et al., 2025 restrictor line paper)
  gam_fits[[ii]] <- 
    fit_mgcv(
      x = cpue_wide,
      common_name = analysis_species[[ii]],
      formula = log_172 ~ s(log_44, bs = "tp", k = 4), 
      family = gaussian(),
      loocv = TRUE)
  
  gam_pred <- dplyr::bind_rows(gam_pred, gam_fits[[ii]]$fit)
  
  gam_boot <- 
    lapply(
      X = boot_cpue_samples, 
      FUN = fit_mgcv, 
      common_name = analysis_species[[ii]],
      formula = log_172 ~ s(log_44, bs = "tp", k = 4), 
      family = gaussian(),
      loocv = FALSE
    ) |>
    lapply(FUN = function(x) x[["fit"]]) %>%
    do.call(rbind, .) |>
    dplyr::mutate(method = "GAM")
  
  loocv_pred <- 
    gam_fits[[ii]]$loocv_results |>
    dplyr::mutate(
      method = "GAM",
      PREDICTED_CPUE_KGKM2_172 = (exp(fit)-1) # Fit
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )
  
  # Kappenman ratio estimator (Kappenman 1992, Wilderbuer et al., 1998, von Szalay and Brown, 2001)
  fpc_out <-
    est_fpc(
      x = cpue_wide,
      common_name = analysis_species[[ii]],
      cpue1_name = "CPUE_KGKM2_44",
      cpue2_name = "CPUE_KGKM2_172",
      method = 4, #Kappenman
      loocv = TRUE,
      nboot = 0
    )

  fpc_fits[[ii]] <- fpc_out$fpc
  
  fpc_boot <- 
    lapply(
      X = boot_cpue_samples, 
      FUN = est_fpc, 
      common_name = analysis_species[[ii]],
      cpue1_name = "CPUE_KGKM2_44",
      cpue2_name = "CPUE_KGKM2_172",
      method = 4, #Kappenman
      loocv = FALSE,
      nboot = 0
    ) |>
    lapply(FUN = function(x) x[["fpc"]]) %>%
    do.call(rbind, .) |>
    dplyr::mutate(method = "K")

  loocv_pred <-
    fpc_out$loocv_results |>
    dplyr::mutate(
      method = "K",
      PREDICTED_CPUE_KGKM2_172 = CPUE_KGKM2_44 * FPC
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )

  # Somerton ratio (Somerton et al., 2002)
  lr_out <- 
    cpue_wide |>
    dplyr::filter(!is.na(CPUE_LOG_RATIO) & !is.infinite(CPUE_LOG_RATIO)) |>
    fit_somerton_ratio(
      common_name = analysis_species[[ii]],
      loocv = TRUE
    )
  
  lr_boot <- 
    lapply(
      X = boot_cpue_samples, 
      FUN = dplyr::filter, 
      !is.na(CPUE_LOG_RATIO) & !is.infinite(CPUE_LOG_RATIO)
    ) |>
    lapply(
      FUN = fit_somerton_ratio, 
      common_name = analysis_species[[ii]],
      loocv = FALSE
    ) |>
    lapply(FUN = function(x) x[["ratio"]]) %>%
    do.call(rbind, .)
  
  lr_boot <- 
    lr_boot |>
    dplyr::select(log_ratio, var, ratio, COMMON_NAME) |>
    dplyr::mutate(method = "LR") |>
    dplyr::bind_rows(lr_boot |>
                       dplyr::select(log_ratio, var, ratio = ratio_bc, COMMON_NAME) |>
                       dplyr::mutate(method = "LR-BC"))
  
  lr_fits[[ii]] <- lr_out$ratio
  
  loocv_pred <- 
    lr_out$loocv_results |>
    dplyr::mutate(
      method = "LR-BC",
      # Bias-corrected ratio * 83-112 CPUE
      PREDICTED_CPUE_KGKM2_172 = CPUE_KGKM2_44*ratio_bc
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )
  
  loocv_pred <- 
    lr_out$loocv_results |>
    dplyr::mutate(
      method = "LR",
      # Bias-corrected ratio * 83-112 CPUE
      PREDICTED_CPUE_KGKM2_172 = CPUE_KGKM2_44*ratio
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )
  
  # CCR mean (O'Leary et a. 2021)
  ccr_out <- 
    est_ccr(
      x = cpue_wide, 
      common_name = analysis_species[[ii]], 
      cpue1_name = "CPUE_KGKM2_172", 
      cpue2_name = "CPUE_KGKM2_44",
      loocv = TRUE,
      bootstrap_ci = FALSE
    )
  
  ccr_boot <- 
    lapply(
      X = boot_cpue_samples, 
      FUN = est_ccr, 
      common_name = analysis_species[[ii]], 
      cpue1_name = "CPUE_KGKM2_172", 
      cpue2_name = "CPUE_KGKM2_44",
      loocv = FALSE,
      bootstrap_ci = FALSE
    ) |>
    lapply(FUN = function(x) x[["ccr"]]) %>%
    do.call(rbind, .) |>
    dplyr::mutate(method = "CCRmean")
  
  ccr_fits[[ii]] <- ccr_out$ccr
  
  loocv_pred <- 
    ccr_out$loocv_results |>
    dplyr::mutate(
      method = "CCRmean",
      PREDICTED_CPUE_KGKM2_172 = CPUE_KGKM2_44*fpc_fit
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )
  
  # CCR beta (Fanning 1984)
  ccr_beta_out <- 
    est_ccr_beta(
    x = cpue_wide, 
    common_name = analysis_species[[ii]], 
    cpue1_name = "CPUE_KGKM2_172", 
    cpue2_name = "CPUE_KGKM2_44",
    loocv = TRUE,
    bootstrap_ci = FALSE
  )
  
  ccr_beta_boot <- 
    lapply(
      X = boot_cpue_samples, 
      FUN = est_ccr_beta, 
      common_name = analysis_species[[ii]], 
      cpue1_name = "CPUE_KGKM2_172", 
      cpue2_name = "CPUE_KGKM2_44",
      loocv = FALSE,
      bootstrap_ci = FALSE
    ) |>
    lapply(FUN = function(x) x[["ccr"]]) %>%
    do.call(rbind, .) |>
    dplyr::mutate(method = "CCRbeta")
  
  ccr_beta_fits[[ii]] <- ccr_beta_out$ccr
  
  loocv_pred <- 
    ccr_beta_out$loocv_results |>
    dplyr::mutate(
      method = "CCRbeta",
      PREDICTED_CPUE_KGKM2_172 = CPUE_KGKM2_44*fpc_fit
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )
  
  
  # Median ratio (O'Leary et al., 2021)
  median_ratio_out <- 
    est_median_ratio(
      x = cpue_wide,
      common_name = analysis_species[[ii]], 
      cpue1_name = "CPUE_KGKM2_172", 
      cpue2_name = "CPUE_KGKM2_44",
      loocv = TRUE, 
      bootstrap_ci = FALSE
    )
  
  median_ratio_boot  <- 
    lapply(
      X = boot_cpue_samples, 
      FUN = est_median_ratio, 
      common_name = analysis_species[[ii]], 
      cpue1_name = "CPUE_KGKM2_172", 
      cpue2_name = "CPUE_KGKM2_44",
      loocv = FALSE,
      bootstrap_ci = FALSE
    ) |>
    lapply(FUN = function(x) x[["ser"]]) %>%
    do.call(rbind, .) |>
    dplyr::mutate(method = "Median")
  
  median_ratio_fits[[ii]] <- median_ratio_out$ser
  
  loocv_pred <- 
    median_ratio_out$loocv_results |>
    dplyr::mutate(
      method = "Median",
      PREDICTED_CPUE_KGKM2_172 = CPUE_KGKM2_44*ser_median
    ) |>
    dplyr::bind_rows(
      loocv_pred
    )
  
  # Bootstrap fits 
  

  
}

# Combine predictions rows and add columns -- NEED TO ADD BOOTSTRAP
fpc_pred <- 
  do.call(rbind, fpc_fits) |>
  dplyr::select(COMMON_NAME, fpc = FPC, fpc_lwr = `Boot_95%_LCI`, fpc_upr = `Boot_95%_UCI`) |>
  dplyr::mutate(method = "K")

ratio_pred_lr <- 
  do.call(rbind, somerton_fits) |>
  dplyr::select(COMMON_NAME, fpc = ratio) |>
  dplyr::mutate(method = "LR")

ratio_pred_lrbc <- 
  do.call(rbind, somerton_fits) |>
  dplyr::select(
    COMMON_NAME, fpc = ratio_bc, 
    fpc_lwr = ratio_lci, 
    fpc_upr = ratio_uci
  ) |>
  dplyr::mutate(method = "LR-BC")

ccr_mean <- 
  do.call(rbind, ccr_fits) |>
  dplyr::select(
    COMMON_NAME, 
    fpc = fpc,
    fpc_lwr = fpc_lwr_ci_boot,
    fpc_upr = fpc_upr_ci_boot
  ) |>
  dplyr::mutate(method = "CCRmean")

ccr_beta <- 
  do.call(rbind, ccr_beta_fits) |>
  dplyr::select(
    COMMON_NAME, 
    fpc = fpc,
    fpc_lwr = fpc_lwr_ci_boot,
    fpc_upr = fpc_upr_ci_boot
  ) |>
  dplyr::mutate(method = "CCRbeta")

median_pred <- 
  do.call(rbind, median_ratio_fits) |>
  dplyr::select(
    COMMON_NAME, 
    fpc = ser_median, 
    fpc_lwr = ser_lwr_ci_boot, 
    fpc_upr = ser_upr_ci_boot
  ) |>
  dplyr::mutate(method = "Median")
  

fpc_est <- bind_rows(
  fpc_pred,
  ratio_pred_lrbc, 
  ratio_pred_lr, 
  ccr_mean,
  ccr_beta,
  median_pred
  )


# GAM Fits

gam_pred$COMMON_NAME <- 
  factor(
    gam_pred$COMMON_NAME, 
    levels = species_codes$COMMON_NAME, 
    labels = species_codes$COMMON_NAME
  )

# Calculate RMSE, total percentage error, and mean percentage error ----
loocv_error_by_haul <- 
  loocv_pred |> 
  dplyr::mutate(
    PREDICTED_WEIGHT_172 = PREDICTED_CPUE_KGKM2_172*AREA_SWEPT_KM2_172,
    PERCENT_ERROR = (PREDICTED_CPUE_KGKM2_172-CPUE_KGKM2_172)/CPUE_KGKM2_172*100,
    ERROR_SQUARED = (PREDICTED_CPUE_KGKM2_172-CPUE_KGKM2_172)^2
  )

loocv_table <- 
  loocv_error_by_haul |>
  dplyr::group_by(
    COMMON_NAME,
    method
  ) |>
  dplyr::summarise(
    SUM_PREDICTED_WEIGHT_172 = sum(PREDICTED_WEIGHT_172, na.rm = TRUE),
    SUM_WEIGHT_172 = sum(WEIGHT_172, na.rm = TRUE),
    MPE = mean(PERCENT_ERROR[!is.infinite(PERCENT_ERROR)], na.rm = TRUE),
    RMSE = sqrt(mean(ERROR_SQUARED, na.rm = TRUE))
    ) |>
  dplyr::mutate(TPE = (SUM_PREDICTED_WEIGHT_172-SUM_WEIGHT_172)/SUM_WEIGHT_172*100)

loocv_long <- 
  loocv_table |>
  dplyr::select(-SUM_WEIGHT_172, -SUM_PREDICTED_WEIGHT_172) |>
  tidyr::pivot_longer(
    cols = c(3, 5)
  ) |>
  dplyr::left_join(
    loocv_table |>
      dplyr::group_by(COMMON_NAME) |>
      dplyr::summarise(RMSE = min(RMSE)) |>
      dplyr::mutate(LOWEST_RMSE = "Yes")
  ) |>
  dplyr::mutate(LOWEST_RMSE = ifelse(is.na(LOWEST_RMSE), "No", "Yes"))

# Save outputs to a file
save(
  fpc_pred,
  ratio_pred,
  gam_pred,
  loocv_pred,
  loocv_error_by_haul,
  loocv_table,
  loocv_long,
  cpue_wide,
  analysis_species,
  gear_codes,
  species_codes,
  file = here::here("analysis", "shelf_slope", "output", "fpc_output.rda")
)

# Setup directories for tables and plots

load(here::here("analysis", "shelf_slope", "output", "fpc_output.rda"))

fpc_plot_dir <- here::here("analysis", "shelf_slope", "plots", "fpc_plots")
sr_plot_dir <- here::here("analysis", "shelf_slope", "plots", "sr_plots")
dir.create(fpc_plot_dir, recursive = TRUE)
dir.create(sr_plot_dir, recursive = TRUE)

# Make tables ----

tab_rmse <- 
  loocv_table |>
  dplyr::select(COMMON_NAME, method, RMSE) |>
  tidyr::pivot_wider(names_from = "method", values_from = "RMSE", values_fn = round)

write.csv(
  tab_rmse,
  file = here::here("analysis", "shelf_slope", "plots", "fpc_plots", "fpc_rmse.csv"),
  row.names = FALSE
)

# Make plots ----

# Mean prediction error - |error| > 100 are set to 100 for comparability
p_err_by_species <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(data = loocv_long,
             mapping = aes(x = method, y = ifelse(abs(value) < 100, value, sign(value) * 100), color = name, size = LOWEST_RMSE), shape = 1) +
  facet_wrap(~COMMON_NAME, ncol = 3) +
  scale_color_colorblind(name = "Value") +
  scale_size_discrete(name = "Lowest RMSE", range = c(1.5, 4)) +
  scale_y_continuous(name = "Percent error (%)") + 
  scale_x_discrete(name = "Method") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title = element_text(size = 8))

png(filename = here::here(fpc_plot_dir, paste0("fpc_err_by_species.png")), width = 169, height = 169, res = 300, units = "mm")
print(p_err_by_species)
dev.off()


# Prediction error for individual species; starting with pollock

p_pred_vs_obs_cpue <- vector(mode = "list", length = length(analysis_species))

for(kk in analysis_species) {
  
  p_pred_vs_obs_cpue <-
    ggplot() +
      geom_point(
        data = dplyr::filter(loocv_error_by_haul, COMMON_NAME == kk),
        mapping = aes(x = PREDICTED_CPUE_KGKM2_172, y = CPUE_KGKM2_172)
      ) + 
      ggtitle(kk) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      facet_wrap(~method) +
      scale_x_log10(name = expression('Predicted CPUE ('*kg%.%km^-2*')')) +
      scale_y_log10(name = expression('Observed CPUE ('*kg%.%km^-2*')')) +
      theme_bw()
  
  p_obs_vs_pcterr <-
    ggplot() +
    geom_point(
      data = dplyr::filter(loocv_error_by_haul, COMMON_NAME == kk),
      mapping = aes(x = CPUE_KGKM2_172, y = ifelse(abs(PERCENT_ERROR) < 100, PERCENT_ERROR, sign(PERCENT_ERROR) * 100 ))
    ) + 
    ggtitle(kk) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~method) +
    scale_x_log10(name = expression('Observed CPUE ('*kg%.%km^-2*')')) +
    scale_y_continuous(name = "Prediction error (%)", expand = c(0,0), limits = c(-105, 105)) +
    theme_bw()

  png(filename = here::here(fpc_plot_dir, paste0("fpc_pred_vs_obs_", gsub("[^[:alnum:]]+", "_", kk), ".png")), width = 169, height = 120, res = 300, units = "mm")
  print(p_pred_vs_obs_cpue)
  dev.off()
  
  png(filename = here::here(fpc_plot_dir, paste0("fpc_cpue_vs_pct_err_", gsub("[^[:alnum:]]+", "_", kk), ".png")), width = 169, height = 120, res = 300, units = "mm")
  print(p_obs_vs_pcterr)
  dev.off()

}

p_mpe_tpe_boxplots <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(data = loocv_long,
             mapping = aes(x = method, y = ifelse(abs(value) < 100, value, sign(value) * 100)),
             outliers = FALSE) +
  geom_point(data = loocv_long,
               mapping = aes(x = method, y = ifelse(abs(value) < 100, value, sign(value) * 100), color = COMMON_NAME), shape = 19, size = rel(3)) +
  scale_y_continuous(name = "Percent error (%)") + 
  scale_color_tableau(name = "Common name", palette = "Tableau 20") +
  facet_wrap(~name, nrow = 2) +
  theme_bw()

png(filename = here::here(fpc_plot_dir, paste0("fpc_tpe_mpe_boxplot.png")), width = 169, height = 120, res = 300, units = "mm")
print(p_mpe_tpe_boxplots)
dev.off()


p_mean_rfp <- ggplot() +
  geom_vline(xintercept = 1, linetype = 3) +
  geom_point(data = fpc_est, mapping = aes(y = COMMON_NAME, x = fpc, color = method),
             size = rel(2.2),
             alpha = 0.7) +
  scale_x_continuous(name = expression('Relative fishing power, '*q[PNE]/q[83-112])) +
  scale_color_tableau() +
  ggtitle("temporary") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "inside",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = element_blank(),
        legend.position.inside = c(0.8, 0.8))

png(filename = here::here(fpc_plot_dir, paste0("fpc_mean_rfp_by_method.png")), width = 120, height = 120, res = 300, units = "mm")
print(p_mean_rfp)
dev.off()

ggplot() +
  geom_point(data = loocv_pred,
             mapping = 
               aes(x = exp(fit)-1, y = ((exp(fit)-1)-CPUE_KGKM2_172)/(CPUE_KGKM2_172)*100
                   ) 
             )+
  geom_abline(intercept = 0, slope = 0, linetype = 2) +
  facet_wrap(~COMMON_NAME) + 
  scale_x_log10(name = expression(CPUE*" ("*kg%.%km^-2*")")) +
  scale_y_continuous(name = "LOOCV relative prediction error (%)") +
  theme_bw()

