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

cpue_1998 <- 
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
    best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
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
  output$fpc$model_name <- aic_table$model_name[aic_table$best]
  
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
          
          sq_errors_median[jj] <- (dat$CPUE_NO_KM2_15[jj] / fpc$ratio - dat$CPUE_NO_KM2_30[jj])^2
          
          sq_errors_mean[jj] <- (dat$CPUE_NO_KM2_15[jj] / fpc$ratio_bc - dat$CPUE_NO_KM2_30[jj])^2
          
          mean_fit[jj] <- dat$COUNT_15[jj] / fpc$ratio_bc * dat$AREA_SWEPT_KM2_30[jj]/dat$AREA_SWEPT_KM2_15[jj]
          
          median_fit[jj] <- dat$COUNT_15[jj] / fpc$ratio * dat$AREA_SWEPT_KM2_30[jj]/dat$AREA_SWEPT_KM2_15[jj]
          
        }
        
      }
      
      results_mean$rmse[ii] <- sqrt(mean(sq_errors_mean))
      results_mean$tpe[ii] <- 100 * (sum(mean_fit)-sum(dat$COUNT_30))/sum(dat$COUNT_30)
      
      results_median$rmse[ii] <- sqrt(mean(sq_errors_median))
      results_median$tpe[ii] <- 100 * (sum(median_fit)-sum(dat$COUNT_30))/sum(dat$COUNT_30)
      
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
    best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table
  )
  
  fpc <-
    data.frame(
      model_name = loocv_table$model_name[loocv_table$best],
      log_ratio = coef(output[['best_model']])['(Intercept)'],
      var = vcov(output[['best_model']])[1,1]
    )
  
  fpc$ratio <- exp(fpc$log_ratio)
  fpc$ratio_lci <- exp(fpc$log_ratio-sqrt(2*fpc$var))
  fpc$ratio_uci <- exp(fpc$log_ratio+sqrt(2*fpc$var))
  
  output$fpc <- fpc
  
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
        
      sq_errors[jj] <- (dat$CPUE_NO_KM2_15[jj] / fpc  - dat$CPUE_NO_KM2_30[jj])^2
      
      fit[jj] <- dat$COUNT_15[jj] / fpc * dat$AREA_SWEPT_KM2_30[jj]/dat$AREA_SWEPT_KM2_15[jj]

      
    }
    
    results$rmse[ii] <- sqrt(mean(sq_errors))
    results$tpe[ii] <- 100 * (sum(fit)-sum(dat$COUNT_30))/sum(dat$COUNT_30)
    
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
    best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table
  )
  
  fpc <-
    data.frame(
      model_name = loocv_table$model_name[loocv_table$best],
      log_ratio = coef(output[['best_model']])['(Intercept)'],
      var = vcov(output[['best_model']])[1,1]
    )
  
  fpc$ratio <- exp(fpc$log_ratio)
  fpc$ratio_lci <- exp(fpc$log_ratio-sqrt(2*fpc$var))
  fpc$ratio_uci <- exp(fpc$log_ratio+sqrt(2*fpc$var))

  output$fpc <- fpc
  
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
    geom_point() +
    geom_smooth(method = 'loess') +
    scale_y_continuous(name = "|Std. residual|") +
    scale_fn(name = x_axis_name) +
    theme_bw()
  
  return(p1)
  
}


# Fitting to Somerton et al. (2002) data from 1998 -------------------------------------------------

fpc_table_1998 <- data.frame()
loocv_table_1998 <- data.frame()

# Snow crab

species_code <- 68580
dat <- cpue_1998[cpue_1998$SPECIES_CODE == species_code, ]

ols_results <- fit_ols(dat)
ols_results$anderson_darling
ols_results$kurtosis
ols_results$cor_test

check_ols_heteroskedasticity(dat = dat, model = ols_results$best_model, predictor = "CPUE_NO_KM2_15", x_axis_name = expression(CPUE[15]*' (#/'*km^2*')'))
check_ols_heteroskedasticity(dat = dat, model = ols_results$best_model, predictor = "CPUE_NO_KM2_30", x_axis_name = expression(CPUE[30]*' (#/'*km^2*')'))

beta_results <- fit_betareg(dat)

binomial_results <- fit_binomial(dat)

fpc_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$fpc,
    ols_results$fpc,
    beta_results$fpc
  ) |>
  dplyr::mutate(common_name = "snow crab") |>
  dplyr::bind_rows(fpc_table_1998)

loocv_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$loocv_table,
    ols_results$loocv_table,
    beta_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = "snow crab") |>
  dplyr::bind_rows(loocv_table_1998)

# par(mfrow = c(2,2))
# plot(beta_results$best_model)
# plot(binomial_results$best_model)

# Tanner crab

species_code <- 68560
dat <- cpue_1998[cpue_1998$SPECIES_CODE == species_code, ]
ols_results <- fit_ols(dat)
ols_results$anderson_darling
ols_results$kurtosis
ols_results$cor_test
check_ols_heteroskedasticity(dat = dat, model = ols_results$best_model, predictor = "CPUE_NO_KM2_15", x_axis_name = expression('CPUE (#/'*km^2*')'))
check_ols_heteroskedasticity(dat = dat, model = ols_results$best_model, predictor = "CPUE_NO_KM2_30", x_axis_name = expression('CPUE (#/'*km^2*')'))

beta_results <- fit_betareg(dat)

binomial_results <- fit_binomial(dat)


fpc_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$fpc,
    ols_results$fpc,
    beta_results$fpc
  ) |>
  dplyr::mutate(common_name = "Tanner crab") |>
  dplyr::bind_rows(fpc_table_1998)

loocv_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$loocv_table,
    ols_results$loocv_table,
    beta_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = "Tanner crab") |>
  dplyr::bind_rows(loocv_table_1998)


par(mfrow = c(2,2))
plot(beta_results$best_model)
plot(binomial_results$best_model)

# Red king crab

species_code <- 69322
dat <- cpue_1998[cpue_1998$SPECIES_CODE == species_code, ]
ols_results <- fit_ols(dat)
ols_results$anderson_darling
ols_results$kurtosis
ols_results$cor_test
ols_results$fpc
ols_results$loocv_table
check_ols_heteroskedasticity(dat = dat, model = ols_results$best_model, predictor = "CPUE_NO_KM2_15", x_axis_name = expression('CPUE (#/'*km^2*')'))
check_ols_heteroskedasticity(dat = dat, model = ols_results$best_model, predictor = "CPUE_NO_KM2_30", x_axis_name = expression('CPUE (#/'*km^2*')'))

beta_results <- fit_betareg(dat)
beta_results$loocv_table

binomial_results <- fit_binomial(dat)
binomial_results$loocv_table

fpc_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$fpc,
    ols_results$fpc,
    beta_results$fpc
  ) |>
  dplyr::mutate(common_name = "red king crab") |>
  dplyr::bind_rows(fpc_table_1998)


loocv_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$loocv_table,
    ols_results$loocv_table,
    beta_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = "red king crab") |>
  dplyr::bind_rows(loocv_table_1998)


par(mfrow = c(2,2))
plot(beta_results$best_model)
plot(binomial_results$best_model)



# Fitting to the full data set ---------------------------------------------------------------------


