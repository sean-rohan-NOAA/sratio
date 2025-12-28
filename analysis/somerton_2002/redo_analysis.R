library(sratio)
library(nortest)

# Replicate Somerton analysis ----

somerton_crab <- 
  readRDS(
    file = here::here("analysis", "somerton_2002", "data", "somerton_crab.rds")
  )

# Load data from Somerton's archived files
somerton_catch <- 
  rbind(
    read.table(file = here::here("analysis", "somerton_2002", "data", "Cb2.txt"),
               header = TRUE, na.strings = ".") |>
      dplyr::mutate(SPECIES_CODE = 68560),
    read.table(file = here::here("analysis", "somerton_2002", "data", "Co2.txt"),
               header = TRUE, na.strings = ".") |>
      dplyr::mutate(SPECIES_CODE = 68580),
    read.table(file = here::here("analysis", "somerton_2002", "data", "RK2.txt"),
               header = TRUE, na.strings = ".") |>
      dplyr::mutate(SPECIES_CODE = 69322)
  ) |>
  dplyr::mutate(CATCH15F = CATCH15T-CATCH15M,
                CATCH30F = CATCH30T-CATCH30M) |>
  dplyr::select(TOW_PAIR = OBS, SPECIES_CODE, CATCH15M, CATCH30M, CATCH15F, CATCH30F) |>
  tidyr::pivot_longer(cols = c("CATCH15M", "CATCH30M", "CATCH15F", "CATCH30F")) |>
  dplyr::mutate(TREATMENT = factor(ifelse(stringr::str_detect(name, "30"), 30, 15)),
                SEX = factor(ifelse(stringr::str_detect(name, "F"), "F", "M"))) |>
  dplyr::select(-name) |>
  tidyr::pivot_wider(values_from = "value", names_from = "TREATMENT", names_prefix = "COUNT_")

somerton_effort <- 
  read.table(file = here::here("analysis", "somerton_2002", "data", "Cb1.txt"),
             header = TRUE, na.strings = ".") |>
  dplyr::select(TOW_PAIR = OBS, EFFORT15, EFFORT30, VESSEL) |>
  tidyr::pivot_longer(cols = c("EFFORT15", "EFFORT30")) |>
  dplyr::mutate(TREATMENT = factor(ifelse(stringr::str_detect(name, "30"), 30, 15))) |>
  dplyr::select(-name) |>
  dplyr::mutate(value = value/100,
                VESSEL = factor(VESSEL)) |>
  tidyr::pivot_wider(values_from = "value", names_from = "TREATMENT", names_prefix = "AREA_SWEPT_KM2_")

cpue <- 
  dplyr::inner_join(somerton_catch, somerton_effort) |>
  dplyr::mutate(
    CPUE_NO_KM2_15 = COUNT_15 / AREA_SWEPT_KM2_15,
    CPUE_NO_KM2_30 = COUNT_30 / AREA_SWEPT_KM2_30,
    LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
    LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
    CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
    CPUE_RATIO = CPUE_NO_KM2_15/CPUE_NO_KM2_30,
    COMBINED_COUNT = COUNT_30 + COUNT_15,
    PROP_15 = COUNT_15/(COUNT_15+COUNT_30),
    EFFORT_RATIO = AREA_SWEPT_KM2_15/AREA_SWEPT_KM2_30,
    common_name = sratio::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0)

# CPUE ratio regression ----- 


make_aic_table <- 
  function(model_list) {
  
  if(is.null(names(model_list))) {
    model_names <- 1:length(model_list)
  } else {
    model_names <- names(model_list)
  }
  
  results <- data.frame(
    model_name = model_names,
    formula = sapply(model_list, function(m) format(formula(m))),
    aic = round(sapply(model_list, AIC), 2),
    k = sapply(model_list, function(m) length(coef(m))),
    stringsAsFactors = FALSE
  )
  
  min_aic <- min(results$aic)
  results$delta_aic <- results$aic - min_aic
  
  candidate_indices <- which(results$delta_aic < 2)
  candidates <- results[candidate_indices, ]
  
  best_row_idx <- candidate_indices[which.min(candidates$k)]
  
  results$best <- FALSE
  results$best[best_row_idx] <- TRUE
  
  return(results[order(results$aic), ])
  
} 

fit_ols <- 
  function(x) {
  
  # Fit models
  ols1 <- 
    lm(
      formula = CPUE_LOG_RATIO ~ 1,
      data = x
    )
  
  # ols2 <- 
  #   lm(
  #     formula = CPUE_LOG_RATIO ~ VESSEL,
  #     data = x
  #   )
  # 
  # ols3 <- 
  #   lm(
  #     formula = CPUE_LOG_RATIO ~ SEX,
  #     data = x
  #   )
  # 
  # ols4 <- 
  #   lm(
  #     formula = CPUE_LOG_RATIO ~ SEX + VESSEL,
  #     data = x
  #   )
  
  model_list <- 
    list(
      ols1 = ols1#,
      # ols2 = ols2,
      # ols3 = ols3,
      # ols4 = ols4
    )
  
  aic_table <- 
    make_aic_table(
      model_list = model_list
  )
  
  loocv_table <- 
    loocv_ols(
      model_list = model_list,
      dat = dat
    )
  
  output <- list(
    best_model = model_list[[aic_table$model_name[aic_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table
  )
  
  output[['cor_test']] <- list(
    CPUE_30 = summary(lm(abs(rstandard(output[['best_model']]))~I(log(x$CPUE_NO_KM2_30)))),
    CPUE_15 = summary(lm(abs(rstandard(output[['best_model']]))~I(log(x$CPUE_NO_KM2_15))))
  )
  
  output[['anderson_darling']] <- nortest::ad.test(rstandard(output[['best_model']]))
  output[['kurtosis']] <- e1071::kurtosis(rstandard(output[['best_model']]), type = 2)
  
  output$fpc <- miller_bias_correct(output[['best_model']])
  
  # Make diagnostic plots for the best-fit model
  return(output)
  
  }

# LOOCV
loocv_ols <- 
  function(model_list, dat) {
    
    results_mean <- results_median <- data.frame(
      model_name = names(model_list),
      formula = sapply(model_list, function(m) format(formula(m))),
      rmse = numeric(length(model_list))
    )
    
    results_mean$method <- "OLS mean"
    results_median$method <- "OLS median"
    
    for(ii in seq_along(model_list)) {
      mod <- model_list[[ii]]
      n_obs <- nrow(dat)
      
      sq_errors_mean <- rep(Inf, n_obs)
      sq_errors_median <- rep(Inf, n_obs)
      
      mean_fit <- rep(Inf, n_obs)
      median_fit <- rep(Inf, n_obs)
      
      
      for(jj in 1:n_obs) {
        
        train_dat <- dat[-jj, , drop = FALSE]
        test_dat  <- dat[jj, , drop = FALSE]
        
        fit_loocv <- update(mod, data = train_dat)
        
        if(length(coef(fit_loocv)) == 1) {
          fpc <- miller_bias_correct(fit_loocv)
          
          sq_errors_median[jj] <- (dat$CPUE_NO_KM2_30[jj] * fpc$ratio - dat$CPUE_NO_KM2_15[jj])^2
          
          sq_errors_mean[jj] <- (dat$CPUE_NO_KM2_30[jj] * fpc$ratio_bc - dat$CPUE_NO_KM2_15[jj])^2
          
          mean_fit[jj] <- dat$COUNT_30[jj] * fpc$ratio_bc * dat$AREA_SWEPT_KM2_15[jj]/dat$AREA_SWEPT_KM2_30[jj]
          
          median_fit[jj] <- dat$COUNT_30[jj] * fpc$ratio * dat$AREA_SWEPT_KM2_15[jj]/dat$AREA_SWEPT_KM2_30[jj]
          
        }
        
      }
      
      results_mean$rmse[ii] <- sqrt(mean(sq_errors_mean))
      results_mean$tpe[ii] <- 100 * (sum(mean_fit)-sum(dat$COUNT_15))/sum(dat$COUNT_15)
      
      results_median$rmse[ii] <- sqrt(mean(sq_errors_median))
      results_median$tpe[ii] <- 100 * (sum(median_fit)-sum(dat$COUNT_15))/sum(dat$COUNT_15)
      
    }
    
    results <- rbind(
      results_mean, results_median
    )
    
    results$best <- results$rmse == min(results$rmse)
    
    return(results)
    
  }

# Function to extract model intercept, variance, and bias-corrected ratio from log-ratio models
miller_bias_correct <- 
  function(mod) {
    
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


fit_betareg <- 
  function(x) {
  
  beta1 <- 
    betareg::betareg(
      formula = PROP_15 ~ 1 | 1,
      link = "logit",
      offset = log(EFFORT_RATIO),
      data = x
    )
  
  beta2 <- 
    betareg::betareg(
      formula = PROP_15 ~ 1 | LOG_CPUE_NO_KM2_15,
      link = "logit",
      offset = log(EFFORT_RATIO),
      data =  x
    )
  
  beta3 <- 
    betareg::betareg(
      formula = PROP_15 ~ 1 | LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
      link = "logit",
      offset = log(EFFORT_RATIO),
      data =  x
    )
  
  model_list <- 
    list(
      beta1 = beta1,
      beta2 = beta2,
      beta3 = beta3
    )
  
  aic_table <- 
    make_aic_table(
      model_list 
    )
  
  loocv_table <-
    loocv_betareg_binomial(
      model = model_list,
      dat = x,
      method = "beta"
    )
  
  output <- list(
    best_model = model_list[[aic_table$model_name[aic_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table
  )
  
  # Residual plots
  
  return(output)
  
  }

loocv_betareg_binomial <- function(model_list, dat, method) {
  
  results <- 
    data.frame(
    model_name = names(model_list),
    formula = sapply(model_list, function(m) format(formula(m))),
    rmse = numeric(length(model_list)),
    method = method
  )
  
  for(ii in seq_along(model_list)) {
    mod <- model_list[[ii]]
    n_obs <- nrow(dat)
    
    sq_errors <- rep(Inf, n_obs)
    fit <- rep(Inf, n_obs)
    
    for(jj in 1:n_obs) {
      
      train_dat <- dat[-jj, , drop = FALSE]
      test_dat  <- dat[jj, , drop = FALSE]
      
      fit_loocv <- update(mod, data = train_dat)
      
      fpc <- exp(coef(fit_loocv)['(Intercept)'])
        
      sq_errors[jj] <- (dat$CPUE_NO_KM2_30[jj] * fpc  - dat$CPUE_NO_KM2_15[jj])^2
      
      fit[jj] <- dat$COUNT_30[jj] * fpc * dat$AREA_SWEPT_KM2_15[jj]/dat$AREA_SWEPT_KM2_30[jj]

      
    }
    
    results$rmse[ii] <- sqrt(mean(sq_errors))
    results$tpe[ii] <- 100 * (sum(fit)-sum(dat$COUNT_15))/sum(dat$COUNT_15)
    
  }
  
  results$best <- results$rmse == min(results$rmse)
  
  return(results)
  
}



fit_binomial <- 
  function(x) {
  
  binomial1 <-
    glm(
      formula = PROP_15 ~ 1,
      family = binomial(link = "logit"),
      weight = COMBINED_COUNT,
      offset = log(EFFORT_RATIO),
      data = x
    )
  
  model_list <- 
    list(
      binomial1 = binomial1
    )
  
  aic_table <- 
    make_aic_table(
      model_list 
    )
  
  loocv_table <- loocv_betareg_binomial(
    model_list = model_list,
    dat = x,
    method = "binomial"
  )
  
  output <- list(
    best_model = model_list[[aic_table$model_name[aic_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table
  )
  
  return(output)
  
}

check_ols_heteroskedasticity <- 
  function(dat, model, predictor, x_axis_name = "Observation", scale = "log10", residual_fn = rstandard) {
  
  if(scale == "log10") {
    scale_fn <- scale_x_log10
  } else {
    scale_fn <- scale_x_continuous
  }
  
  residual_df <-
    data.frame(
      Observed = dat[[predictor]],
      Residual = abs(residual_fn(model))
    )
  
  p1 <- 
    ggplot(
      data = residual_df,
      mapping = aes(x = Observed, y = Residual)) +
    # geom_hline(yintercept = 0, linetype = 2) +
    geom_point() +
    geom_smooth(method = 'loess') +
    scale_y_continuous(name = "|Std. residual|") +
    scale_fn(name = x_axis_name) +
    theme_bw()
  
  return(p1)
  
}



species_code <- 68580
dat <- cpue[cpue$SPECIES_CODE == species_code, ]
somerton_results <- fit_ols(dat)
somerton_results$anderson_darling
somerton_results$kurtosis
somerton_results$cor_test
somerton_results$fpc
somerton_results$loocv_table
check_ols_heteroskedasticity(dat = dat, model = somerton_results$best_model, predictor = "CPUE_NO_KM2_15", x_axis_name = expression('CPUE (#/'*km^2*')'))
check_ols_heteroskedasticity(dat = dat, model = somerton_results$best_model, predictor = "CPUE_NO_KM2_30", x_axis_name = expression('CPUE (#/'*km^2*')'))
beta_results <- fit_betareg(dat)
beta_results$loocv_table
binomial_results <- fit_binomial(dat)

par(mfrow = c(2,2))
plot(beta_results$best_model)
plot(binomial_results$best_model)

hist(resid(beta_results$best_model))

e1071::kurtosis(resid(beta_results$best_model))
e1071::kurtosis(rstandard(somerton_results$best_model))


species_code <- 68560
dat <- cpue[cpue$SPECIES_CODE == species_code, ]
somerton_results <- fit_ols(dat)
somerton_results$anderson_darling
somerton_results$kurtosis
somerton_results$cor_test
somerton_results$fpc
somerton_results$loocv_table
check_ols_heteroskedasticity(dat = dat, model = somerton_results$best_model, predictor = "CPUE_NO_KM2_15", x_axis_name = expression('CPUE (#/'*km^2*')'))
check_ols_heteroskedasticity(dat = dat, model = somerton_results$best_model, predictor = "CPUE_NO_KM2_30", x_axis_name = expression('CPUE (#/'*km^2*')'))
beta_results <- fit_betareg(dat)
binomial_results <- fit_binomial(dat)

par(mfrow = c(2,2))
plot(beta_results$best_model)
plot(binomial_results$best_model)



species_code <- 69322
dat <- cpue[cpue$SPECIES_CODE == species_code, ]
somerton_results <- fit_ols(dat)
somerton_results$anderson_darling
somerton_results$kurtosis
somerton_results$cor_test
somerton_results$fpc
somerton_results$loocv_table
check_ols_heteroskedasticity(dat = dat, model = somerton_results$best_model, predictor = "CPUE_NO_KM2_15", x_axis_name = expression('CPUE (#/'*km^2*')'))
check_ols_heteroskedasticity(dat = dat, model = somerton_results$best_model, predictor = "CPUE_NO_KM2_30", x_axis_name = expression('CPUE (#/'*km^2*')'))
beta_results <- fit_betareg(dat)
binomial_results <- fit_binomial(dat)


par(mfrow = c(2,2))
plot(beta_results$best_model)
plot(binomial_results$best_model)


# Make lines
one_to_one <- 
  data.frame(SPECIES_CODE = c(69322, 68580, 68560),
             slope = 1,
             intercept = 0)

# Values from Table 2 in Somerton et al. (2002)
somerton_reported <-
  data.frame(
    SPECIES_CODE = c(69322, 68580, 68560),
    common_name = sratio::species_code_label(x = c(69322, 68580, 68560), type = "common_name"),
    ratio_bc = c(1.244, 1.784, 1.681),
    var = c(0.401, 0.626, 0.583)
  )

somerton_reported$ratio <- exp(log(somerton_reported$ratio_bc) - 0.5 * somerton_reported$var)

pred_df <-
  rbind(
    data.frame(
      CPUE_NO_KM2_30 = 
        exp(seq(min(log(cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560])), max(log(cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560])), by = 0.1)),
      SPECIES_CODE = 68560
    ),
    data.frame(
      CPUE_NO_KM2_30 = 
        exp(seq(min(log(cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580])), max(log(cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580])), by = 0.1)),
      SPECIES_CODE = 68580
    ),
    data.frame(
      CPUE_NO_KM2_30 = 
        exp(seq(min(log(cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322])), max(log(cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322])), by = 0.1)),
      SPECIES_CODE = 69322
    )
  ) |>
  dplyr::inner_join(somerton_reported) |>
  dplyr::mutate(fit_bc = CPUE_NO_KM2_30*ratio_bc,
                fit = CPUE_NO_KM2_30*ratio)

p2 <- 
  ggplot() +
  geom_abline(
    data = one_to_one,
    mapping = aes(color = "1:1 line",
                  linetype = "1:1 line",
                  slope = slope,
                  intercept = intercept)
  ) +
  geom_path(
    data = pred_df,
    mapping = 
      aes(x = CPUE_NO_KM2_30, y = fit_bc, 
          color = "Somerton BC", 
          linetype = "Somerton BC")
  ) +
  geom_path(
    data = pred_df,
    mapping = 
      aes(x = CPUE_NO_KM2_30, y = fit, 
          color = "Somerton raw", 
          linetype = "Somerton raw")
  ) +
  geom_line(
    data = zeroint_fit,
    mapping = aes(x = exp(LOG_CPUE_NO_KM2_30),
                  y = fit,
                  linetype = "lm raw",
                  color = "lm raw")
  ) +
  geom_line(
    data = zeroint_fit,
    mapping = aes(x = exp(LOG_CPUE_NO_KM2_30),
                  y = fit_bc,
                  linetype = "lm BC",
                  color = "lm BC")
  ) +
  geom_point(
    data = cpue,
    mapping = aes(x = CPUE_NO_KM2_30, y = CPUE_NO_KM2_15),
    size = rel(0.6)
  ) +
  scale_color_manual(name = NULL, 
                     values = 
                       c("1:1 line" = "black",
                         "Somerton BC" = "red",
                         "Somerton raw" = "salmon",
                         "lm BC" = "darkgreen",
                         "lm raw" = "green"
                       )
  ) +
  scale_linetype_manual(name = NULL, 
                        values = 
                          c("1:1 line" = 2,
                            "Somerton BC" = 1,
                            "Somerton raw" = 1,
                            "lm BC" = 1,
                            "lm raw" = 1
                          )
  ) +
  scale_x_log10(name = expression(CPUE[30]~('#'%.%km^-2))) +
  scale_y_log10(name = expression(CPUE[15]~('#'%.%km^-2))) +
  facet_wrap(~common_name, scales = "free") +
  theme_bw() +
  theme(legend.title = element_blank())

print(p2)

png(filename = here::here("analysis", "somerton_2002", "plots", "compare_somerton.png"),
    width = 8, height = 3, units = "in", res = 300)
print(p2)
dev.off()

write.csv(
  ratios, 
  file = here::here("analysis", "somerton_2002", "plots", "ratio_estimates.csv"), 
  row.names = FALSE
)
