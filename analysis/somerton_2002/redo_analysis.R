library(sratio)
library(nortest)
library(xlsx)

# Prep data -- move this to a separate file

# Somerton's archived data

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

saveRDS(object = cpue_1998, file = here::here("analysis", "somerton_2002", "data", "cpue_1998.rds"))


# 1995 + 2021-2024 experiments

cpue_other <- 
  sratio::data_1530$size |>
  dplyr::filter(SPECIES_CODE %in% c(68580, 68560, 69322)) |>
  dplyr::mutate(
    SEX = ifelse(SEX == 1, "M", "F"),
    MATCHUP = paste0("O", MATCHUP)
                ) |>
  dplyr::group_by(VESSEL, CRUISE, HAUL, SPECIES_CODE, SEX) |>
  dplyr::summarise(COUNT = sum(SAMPLING_FACTOR)) |>
  dplyr::inner_join(
    dplyr::select(
      data_1530$haul,
      VESSEL, CRUISE, HAUL, AREA_SWEPT_KM2, TREATMENT, YEAR, TOW_PAIR = MATCHUP
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(YEAR != 1998) |>
  dplyr::select(-VESSEL, -CRUISE, -HAUL) |>
  tidyr::pivot_wider(values_from = c("COUNT", "AREA_SWEPT_KM2"), names_from = TREATMENT, values_fill = 0) |>
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

saveRDS(object = cpue_other, file = here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))
  

# Load data ----

cpue_other <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))
cpue_1998 <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_1998.rds"))


# CPUE ratio regression ----- 

make_aic_table <- 
  function(model_list) {
    
    results <- data.frame(
      model_name = names(model_list) %||% seq_along(model_list),
      formula    = sapply(model_list, function(m) paste(format(formula(m)), collapse = "")),
      aic        = round(sapply(model_list, AIC), 2),
      k          = sapply(model_list, function(m) attr(logLik(m), "df")),
      stringsAsFactors = FALSE
    )
    
    results$delta_aic <- results$aic - min(results$aic)
    
    candidates <- results[results$delta_aic < 2, ]
    best_name <- candidates$model_name[which.min(candidates$k)]
    results$best <- results$model_name == best_name
    
    return(results[order(results$aic), ])
  
} 



fit_ols <- 
  function(x, bootstrap_samples = NULL) {

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
  
  bootstrap_results <- NA
  
  if(!is.null(bootstrap_samples)) {
    bootstrap_results <- boot_ols(
      model_list = model_list,
      boot_samples_list = bootstrap_samples
    )
  }
  
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
    models = model_list,
    best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table,
    bootstrap_results = bootstrap_results
  )
  
  output[['cor_test']] <- list(
    CPUE_30 = summary(lm(abs(rstandard(output[['best_model']]))~I(log(x$CPUE_NO_KM2_30)))),
    CPUE_15 = summary(lm(abs(rstandard(output[['best_model']]))~I(log(x$CPUE_NO_KM2_15))))
  )
  
  output[['anderson_darling']] <- nortest::ad.test(rstandard(output[['best_model']]))
  output[['kurtosis']] <- e1071::kurtosis(rstandard(output[['best_model']]), type = 2)
  
  fpc <- miller_bias_correct(output[['best_model']])
  fpc$model_name <- aic_table$model_name[aic_table$best]
  
  fpc <- 
    dplyr::bind_rows(
      fpc |>
        dplyr::select(-ratio, -ratio_lci, -ratio_uci) |>
        dplyr::rename(ratio = ratio_bc, ratio_lci = ratio_bc_lci, ratio_uci = ratio_bc_uci) |>
        dplyr::mutate(method = "OLS mean"),
      fpc |>
        dplyr::select(-ratio_bc, -ratio_bc_lci, -ratio_bc_uci) |>
        dplyr::mutate(method = "OLS median")
    )
  
  output$fpc <- fpc
  
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
          
          mean_fit[jj] <- dat$CPUE_NO_KM2_15[jj] / fpc$ratio_bc * dat$AREA_SWEPT_KM2_30[jj]
          
          median_fit[jj] <- dat$CPUE_NO_KM2_15[jj] / fpc$ratio * dat$AREA_SWEPT_KM2_30[jj]
          
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


boot_ols <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      mod <- model_list[[ii]]
      
      fits <- data.frame()
      
      for(jj in 1:length(boot_samples_list)) {
        
        boot_fit <- update(mod, data = boot_samples_list[[jj]])
        
        boot_fpc <- miller_bias_correct(mod = boot_fit)
        
        boot_fpc$model_name <- names(model_list)[ii]
        
        fits <- dplyr::bind_rows(
          fits,
          boot_fpc
        )
        
      }
      
      model_fits[[ii]] <- fits
      
      bootstrap_median <-
        data.frame(
          model_name = names(model_list)[ii],
          method = "median",
          
          ratio     = quantile(fits$ratio, 0.5),
          ratio_lci        = quantile(fits$ratio, 0.025),
          ratio_uci        = quantile(fits$ratio, 0.975)
        )
      
      bootstrap_mean <-
        data.frame(
          model_name = names(model_list)[ii],
          method = "mean",
          
          ratio  = quantile(fits$ratio_bc, 0.5),
          ratio_lci     = quantile(fits$ratio_bc, 0.025),
          ratio_uci     = quantile(fits$ratio_bc, 0.975)
        )
      
      bootstrap_results <-
        dplyr::bind_rows(
          bootstrap_median,
          bootstrap_mean,
          bootstrap_results
        )
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
  }


# Function to extract model intercept, variance, and bias-corrected ratio from log-ratio models
miller_bias_correct <- 
  function(mod) {
    
    if(all(names(coef(mod)) == "(Intercept)")) {
      
      log_ratio <- coef(mod)["(Intercept)"]
      var <- summary(mod)$sigma^2
      se <- summary(mod)$coefficients[, 2]
      ratio <- exp(log_ratio)
      ratio_lci <- exp(log_ratio - 2 * se)
      ratio_uci <- exp(log_ratio + 2 * se)
      ratio_bc <- exp(log_ratio + 0.5 * var)
      ratio_bc_lci <- exp(log_ratio - 2 * se + 0.5 * var)
      ratio_bc_uci <- exp(log_ratio + 2 * se + 0.5 * var)
      
      output <- data.frame(
        log_ratio = unname(log_ratio),
        var = unname(var),
        ratio = unname(ratio),
        ratio_lci = ratio_lci,
        ratio_uci = ratio_uci,
        ratio_bc = unname(ratio_bc),
        ratio_bc_lci = ratio_bc_lci,
        ratio_bc_uci = ratio_bc_uci
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
  function(x, bootstrap_samples = NULL) {
      
    beta1 <- 
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        family = beta_family(link = "logit"),
        offset = log(EFFORT_RATIO),
        data = x
      )
  
    beta2 <- 
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        dispformula = ~ LOG_CPUE_NO_KM2_15,
        family = beta_family(link = "logit"),
        offset = log(EFFORT_RATIO),
        data = x
      )
    
    beta3 <- 
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        dispformula = ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
        family = beta_family(link = "logit"),
        offset = log(EFFORT_RATIO),
        data = x
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
    loocv_betareg(
      model = model_list,
      dat = x,
      method = "beta"
    )
  
  bootstrap_results <- NA
  
  if(!is.null(bootstrap_samples)) {
    bootstrap_results <- boot_betareg_binomial(
      model_list = model_list,
      boot_samples_list = bootstrap_samples
    )
  }
  
  output <- list(
    models = model_list,
    best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
    aic_table = aic_table,
    loocv_table = loocv_table,
    bootstrap_results = bootstrap_results
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



loocv_betareg <- 
  function(model_list, dat, method) {
    
    results <- 
      data.frame(
        model_name = names(model_list),
        formula = sapply(model_list, function(m) format(formula(m))),
        rmse = numeric(length(model_list)),
        method = method
      )
    
    for(ii in seq_along(model_list)) {
      
      cat("loocv ", names(model_list)[ii], "\n")
      mod <- model_list[[ii]]
      n_obs <- nrow(dat)
      
      sq_errors <- rep(Inf, n_obs)
      fit_r <- rep(Inf, n_obs)
      
      for(jj in 1:n_obs) {  
        
        train_dat <- dat[-jj, , drop = FALSE]
        test_dat  <- dat[jj, , drop = FALSE]
        
        fit_loocv <- update(mod, data = train_dat)
        
        # Predict proportion
        fit_r[jj] <- exp(fixef(fit_loocv)$cond['(Intercept)'])
        
      }
      
      sq_errors <- (dat$CPUE_NO_KM2_15 /fit_r  - dat$CPUE_NO_KM2_30)^2
      
      fit <- dat$CPUE_NO_KM2_15 / fit_r * dat$AREA_SWEPT_KM2_30
      
      results$rmse[ii] <- sqrt(mean(sq_errors))
      results$tpe[ii] <- 100 * (sum(fit)-sum(dat$COUNT_30))/sum(dat$COUNT_30)
      
    }
    
    results$best <- results$rmse == min(results$rmse)
    
    return(results)
    
  }

boot_betareg <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      cat(names(model_list)[ii], "\n")
      
      mod <- model_list[[ii]]
      
      fits <- data.frame()
      
      for(jj in 1:length(boot_samples_list)) {
        
        if(jj%%100 == 0) {
          cat(jj, "\n")
        }
        
        boot_fit <- update(mod, data = boot_samples_list[[jj]])
        
        coef <- 
          fixef(boot_fit) |> 
          as.list() |>
          do.call(what = rbind) |>
          as.data.frame()
        
        type <- rownames(coef)
        rownames(coef) <- NULL
        coef$type <- type
        coef$model_name <- names(model_list)[ii]
        
        fits <- dplyr::bind_rows(
          fits,
          coef
        )
        
      }
      
      value_cols <- names(fits)[!(names(fits) %in% c("type", "model_name"))]
      
      ci <- 
        fits %>%
        dplyr::group_by(
          dplyr::across(
            dplyr::all_of(c("type", "model_name")))
        ) %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(value_cols),
            list(
              median = ~median(.x, na.rm = TRUE),
              lci = ~quantile(.x, probs = 0.025),
              uci = ~quantile(.x, probs = 0.975)
            ),
            .names = "{.col}_{.fn}"
          ),
          .groups = "drop"
        )
      
      model_fits[[ii]] <- fits
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
  }


fit_binomial <- 
  function(x, bootstrap_samples = NULL) {
    
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
    
    loocv_table <- 
      loocv_betareg_binomial(
        model_list = model_list,
        dat = x,
        method = "binomial"
      )
    
    bootstrap_results <- NA
    
    if(!is.null(bootstrap_samples)) {
      bootstrap_results <- boot_betareg_binomial(
        model_list = model_list,
        boot_samples_list = bootstrap_samples
      )
    }
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      bootstrap_results = bootstrap_results
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


loocv_binomial <- 
  function(model_list, dat, method) {
  
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
      
      fit[jj] <- dat$CPUE_NO_KM2_15[jj] / fpc * dat$AREA_SWEPT_KM2_30[jj]

      
    }
    
    results$rmse[ii] <- sqrt(mean(sq_errors))
    results$tpe[ii] <- 100 * (sum(fit)-sum(dat$COUNT_30))/sum(dat$COUNT_30)
    
  }
  
  results$best <- results$rmse == min(results$rmse)
  
  return(results)
  
}

boot_binomial <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      mod <- model_list[[ii]]
      
      fits <- data.frame()
      
      for(jj in 1:length(boot_samples_list)) {
        
        boot_fit <- update(mod, data = boot_samples_list[[jj]])
        
        boot_fpc <- as.data.frame(t(coef(boot_fit)))
        boot_fpc$ratio <- exp(boot_fpc[['(Intercept)']])
        
        boot_fpc$model_name <- names(model_list)[ii]
        
        fits <- dplyr::bind_rows(
          fits,
          boot_fpc
        )
        
      }
      
      model_fits[[ii]] <- fits
      
      bootstrap_results <- 
        data.frame(
          model_name = names(model_list)[ii],
          
          ratio = quantile(fits$ratio, 0.5),
          ratio_lci = quantile(fits$ratio, 0.025),
          ratio_uci = quantile(fits$ratio, 0.975)
        ) |>
        dplyr::bind_rows(
          bootstrap_results
        )
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
  }



fit_pois_nb <- 
  function(x, bootstrap_samples = NULL) {
    
    poisson1 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15,
        family = poisson(link = "log"),
        offset = log(1/EFFORT_RATIO),
        data = x
      )
    
    poisson2 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15 + I(COUNT_15^2),
        family = poisson(link = "log"),
        offset = log(1/EFFORT_RATIO),
        data = x
      )
    
    nbin1 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15,
        family = glmmTMB::nbinom1(link = "log"),
        offset = log(1/EFFORT_RATIO),
        data = x
      )
    
    nbin2 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15,
        family = glmmTMB::nbinom1(link = "log"),
        offset = log(1/EFFORT_RATIO),
        disp = ~ CPUE_NO_KM2_15,
        data = x
      )
    
    nbin3 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15,
        family = glmmTMB::nbinom1(link = "log"),
        offset = log(1/EFFORT_RATIO),
        disp = ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
        data = x
      )
    
    nbin4 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15 + I(COUNT_15^2),
        family = glmmTMB::nbinom1(link = "log"),
        offset = log(1/EFFORT_RATIO),
        data = x
      )
    
    nbin5 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15 + I(COUNT_15^2),
        family = glmmTMB::nbinom1(link = "log"),
        offset = log(1/EFFORT_RATIO),
        disp = ~ LOG_CPUE_NO_KM2_15,
        data = x
      )
    
    nbin6 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ COUNT_15 + I(COUNT_15^2),
        family = glmmTMB::nbinom1(link = "log"),
        offset = log(1/EFFORT_RATIO),
        disp = ~  LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
        data = x
      )
    
    model_list <- 
      list(
        poisson1 = poisson1,
        poisson2 = poisson2,
        nbin1 = nbin1,
        nbin2 = nbin2,
        nbin3 = nbin3,
        nbin4 = nbin4,
        nbin5 = nbin5,
        nbin6 = nbin6
      )
    
    aic_table <- 
      dplyr::bind_rows(
        make_aic_table(
          model_list[grepl(pattern = "poisson", x = names(model_list))]
        ),
        make_aic_table(
          model_list[grepl(pattern = "nbin", x = names(model_list))]
        )
      )
    
    loocv_table <- 
      loocv_pois_nb(
        model_list = model_list,
        dat = x,
        method = ""
      )
    
    bootstrap_results <- NA
    
    if(!is.null(bootstrap_samples)) {
      
      cat("Starting bootstrap fits\n")
      
      bootstrap_results <- boot_pois_nb(
        model_list = model_list,
        boot_samples_list = bootstrap_samples
      )
    }
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      bootstrap_results = bootstrap_results,
      fpc = NA
    )
    
    return(output)
    
  }



loocv_pois_nb <- 
  function(model_list, dat, method) {
    
    results <- 
      data.frame(
        model_name = names(model_list),
        formula = sapply(model_list, function(m) format(formula(m))),
        rmse = numeric(length(model_list)),
        method = method
      )
    
    for(ii in seq_along(model_list)) {
      
      cat(names(model_list)[ii], "\n")
      mod <- model_list[[ii]]
      n_obs <- nrow(dat)
      
      sq_errors <- rep(Inf, n_obs)
      fit <- rep(Inf, n_obs)
      
      for(jj in 1:n_obs) {  
        
        train_dat <- dat[-jj, , drop = FALSE]
        test_dat  <- dat[jj, , drop = FALSE]
        
        fit_loocv <- update(mod, data = train_dat)
        
        # Predict counts
        fit[jj] <- predict(fit_loocv, newdata = test_dat, type = "response")
        
        obs <- test_dat$CPUE_NO_KM2_30
        
        sq_errors[jj] <- (fit[jj]/test_dat$AREA_SWEPT_KM2_30 - obs)^2
        
      }
      
      results$rmse[ii] <- sqrt(mean(sq_errors))
      results$tpe[ii] <- 100 * (sum(fit)-sum(dat$COUNT_30))/sum(dat$COUNT_30)
      
    }
    
    results$best <- results$rmse == min(results$rmse)
    
    return(results)
    
  }

boot_pois_nb <- 
  function(model_list, boot_samples_list) {
    
    model_fits <- vector(mode = "list", length = length(model_list))
    names(model_fits) <- names(model_list)
    
    bootstrap_results <- data.frame()
    
    for(ii in seq_along(model_list)) {
      
      cat(names(model_list)[ii], "\n")
      
      mod <- model_list[[ii]]
      
      fits <- data.frame()
      
      for(jj in 1:length(boot_samples_list)) {
        
        if(jj%%100 == 0) {
          cat(jj, "\n")
        }
        
        boot_fit <- update(mod, data = boot_samples_list[[jj]])
        
        coef <- 
          fixef(boot_fit) |> 
          as.list() |>
          do.call(what = rbind) |>
          as.data.frame()
        
        type <- rownames(coef)
        rownames(coef) <- NULL
        coef$type <- type
        coef$model_name <- names(model_list)[ii]
        
        fits <- dplyr::bind_rows(
          fits,
          coef
        )
        
      }
      
      value_cols <- names(fits)[!(names(fits) %in% c("type", "model_name"))]
      
      ci <- 
        fits %>%
        dplyr::group_by(
          dplyr::across(
            dplyr::all_of(c("type", "model_name")))
          ) %>%
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(value_cols),
            list(
              median = ~median(.x, na.rm = TRUE),
              lci = ~quantile(.x, probs = 0.025),
              uci = ~quantile(.x, probs = 0.975)
              ),
            .names = "{.col}_{.fn}"
          ),
          .groups = "drop"
        )
      
      model_fits[[ii]] <- fits
      
    }
    
    output <- 
      list(
        ci = bootstrap_results,
        fits = model_fits
      )
    
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


# Survey-level validation using a different data set
indpendent_fpc_validation <- 
  function(fpc, validation_dat) {
    
    fpc[, c("rmse", "rmse_lci", "rmse_uci", "tpe", "tpe_lci", "tpe_uci")] <- NA
    
    for(ii in 1:nrow(fpc)) {
      
      fit <- validation_dat$CPUE_NO_KM2_15 / fpc$ratio[ii] * validation_dat$AREA_SWEPT_KM2_30
      fit_lci <- validation_dat$CPUE_NO_KM2_15 / fpc$ratio_lci[ii] * validation_dat$AREA_SWEPT_KM2_30
      fit_uci <- validation_dat$CPUE_NO_KM2_15 / fpc$ratio_uci[ii] * validation_dat$AREA_SWEPT_KM2_30

      fpc$rmse[ii] <- sqrt(mean((validation_dat$CPUE_NO_KM2_15 / fpc$ratio[ii] - validation_dat$CPUE_NO_KM2_30)^2))
      fpc$rmse_lci[ii] <- sqrt(mean((validation_dat$CPUE_NO_KM2_15 / fpc$ratio_lci[ii] - validation_dat$CPUE_NO_KM2_30)^2))
      fpc$rmse_uci[ii] <- sqrt(mean((validation_dat$CPUE_NO_KM2_15 / fpc$ratio_uci[ii] - validation_dat$CPUE_NO_KM2_30)^2))

      fpc$tpe[ii] <- 100 * (sum(fit) - sum(validation_dat$COUNT_30))/sum(validation_dat$COUNT_30)
      fpc$tpe_lci[ii] <- 100 * (sum(fit_lci) - sum(validation_dat$COUNT_30))/sum(validation_dat$COUNT_30)
      fpc$tpe_uci[ii] <- 100 * (sum(fit_uci) - sum(validation_dat$COUNT_30))/sum(validation_dat$COUNT_30)
      
    }
    
    return(fpc)
    
  }



twofold_cv <- 
  function(model_list, validation_dat) {
    
    model_list <- beta_results$models
    
    validation_dat <- dat_oos
    
    for(ii in model_list) {
      
      mod <- model_list[[ii]]
    
      # OLS case
      if(is(mod, "lm")) {
        predictions <- 
          predict(
            mod, 
            type = "response", 
            newdata = validation_dat, 
            se.fit = TRUE
          )
        
        validation_dat$fit <- validation_dat$CPUE_NO_KM2_30 * predictions$fit
        validation_dat$fit_lci <- validation_dat$CPUE_NO_KM2_30 * (predictions$fit - 2*predictions$se.fit)
        validation_dat$fit_uci <- validation_dat$CPUE_NO_KM2_30 * (predictions$fit + 2*predictions$se.fit)
        
      }
      
      if(is(mod, "glm")) {
        
        predictions <- 
          predict(
            object = mod, 
            type = "precision", 
            newdata = validation_dat
          )
        
        predictions <- 
          predict(
            mod, 
            type = "terms", 
            newdata = validation_dat
          )
        
        validation_dat$fit <- validation_dat$CPUE_NO_KM2_30 * predictions$fit
        validation_dat$fit_lci <- validation_dat$CPUE_NO_KM2_30 * (predictions$fit - 2*predictions$se.fit)
        validation_dat$fit_uci <- validation_dat$CPUE_NO_KM2_30 * (predictions$fit + 2*predictions$se.fit)
        
        
      }
      
      if(is(mod, "glmmTMB")) {
        
        
      }

      
      validation_dat$fit <- predictions$fit
      
    }
    
  }



draw_bootstrap_samples <- 
  function(x, seed = NULL, grouping_var, n_draws = 1000, replace = TRUE) {
    
    grouping_var_values <- unique(x[[grouping_var]])
    n_obs <- length(grouping_var_values)
    
    draws <- vector(mode = "list", length = n_draws)
    
    if(is.numeric(seed)) {
      set.seed(seed)
    }
    
    for(ii in 1:n_draws) {
      
      samp <- data.frame(var = sample(x = grouping_var_values, replace = replace, size = n_obs))
      
      draws[[ii]] <-
        dplyr::inner_join(
          x, 
          samp, 
          by = setNames("var", grouping_var), 
          relationship = "many-to-many"
      ) |>
        dplyr::mutate(draw = ii)
    }
    
    return(draws)
    
  }


# Fitting Somerton et al. (2002) 1998 data ---------------------------------------------------------

fpc_table_1998 <- data.frame()
loocv_table_1998 <- data.frame()
oos_table <- data.frame()

# Snow crab

# Set species
species_code <- 68580
common_name <- "snow crab"
seed <- 1337
n_draws <- 100

# Setup data
dat <- cpue_1998[cpue_1998$SPECIES_CODE == species_code, ]
dat_bootstrap <- 
  draw_bootstrap_samples(
    x = dat, 
    seed = seed, 
    grouping_var = "TOW_PAIR", 
    n_draws = n_draws, 
    replace = TRUE
  )
dat_oos <- cpue_other[cpue_other$SPECIES_CODE == species_code, ]

# Fit models
# ols_results <- fit_ols(x = dat, bootstrap_samples = dat_bootstrap)
# beta_results <- fit_betareg(x = dat, bootstrap_samples = dat_bootstrap)
# binomial_results <- fit_binomial(x = dat, bootstrap_samples = dat_bootstrap)
# pois_nb_results <- fit_pois_nb(x = dat, bootstrap_samples = dat_bootstrap)

ols_results <- fit_ols(x = dat)
beta_results <- fit_betareg(x = dat)
binomial_results <- fit_binomial(x = dat)
pois_nb_results <- fit_pois_nb(x = dat)

ols_results$anderson_darling
ols_results$kurtosis
ols_results$cor_test
ols_results$bootstrap_results$ci
ols_results$fpc


p1 <- 
  check_ols_heteroskedasticity(
  dat = dat, 
  model = ols_results$best_model, 
  predictor = "CPUE_NO_KM2_15", 
  x_axis_name = expression(CPUE[15]*' (#/'*km^2*')')
  )

p2 <- 
  check_ols_heteroskedasticity(
  dat = dat, 
  model = ols_results$best_model, 
  predictor = "CPUE_NO_KM2_30", 
  x_axis_name = expression(CPUE[30]*' (#/'*km^2*')')
)


cowplot::plot_grid(
  p1, p2
)



fpc_1998 <- 
  dplyr::bind_rows(
  binomial_results$bootstrap_results$ci,
  ols_results$bootstrap_results$ci,
  beta_results$bootstrap_results$ci
)

fpc_table_1998 <- 
  fpc_1998 |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(fpc_table_1998)

loocv_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$loocv_table,
    ols_results$loocv_table,
    beta_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(loocv_table_1998)

oos_table <- 
  indpendent_fpc_validation(
    fpc = fpc_1998,
    validation_dat = dat_oos
  ) |> 
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(oos_table)
  

# par(mfrow = c(2,2))
# plot(beta_results$best_model)
# plot(binomial_results$best_model)

# Tanner crab

species_code <- 68560
common_name <- "Tanner crab"

# Setup data
dat <- cpue_1998[cpue_1998$SPECIES_CODE == species_code, ]
dat_bootstrap <- 
  draw_bootstrap_samples(
    x = dat, 
    seed = seed, 
    grouping_var = "TOW_PAIR", 
    n_draws = 1000, 
    replace = TRUE
  )
dat_oos <- cpue_other[cpue_other$SPECIES_CODE == species_code, ]

# Fit models
ols_results <- fit_ols(x = dat, bootstrap_samples = dat_bootstrap)
beta_results <- fit_betareg(x = dat, bootstrap_samples = dat_bootstrap)
binomial_results <- fit_binomial(x = dat, bootstrap_samples = dat_bootstrap)

ols_results$anderson_darling
ols_results$kurtosis
ols_results$cor_test
ols_results$bootstrap_results$ci
ols_results$fpc


p1 <- 
  check_ols_heteroskedasticity(
    dat = dat, 
    model = ols_results$best_model, 
    predictor = "CPUE_NO_KM2_15", 
    x_axis_name = expression(CPUE[15]*' (#/'*km^2*')')
  )

p2 <- 
  check_ols_heteroskedasticity(
    dat = dat, 
    model = ols_results$best_model, 
    predictor = "CPUE_NO_KM2_30", 
    x_axis_name = expression(CPUE[30]*' (#/'*km^2*')')
  )


cowplot::plot_grid(
  p1, p2
)



fpc_1998 <- 
  dplyr::bind_rows(
    binomial_results$bootstrap_results$ci,
    ols_results$bootstrap_results$ci,
    beta_results$bootstrap_results$ci
  )

fpc_table_1998 <- 
  fpc_1998 |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(fpc_table_1998)

loocv_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$loocv_table,
    ols_results$loocv_table,
    beta_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(loocv_table_1998)

oos_table <- 
  indpendent_fpc_validation(
    fpc = fpc_1998,
    validation_dat = dat_oos
  ) |> 
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(oos_table)


# par(mfrow = c(2,2))
# plot(beta_results$best_model)
# plot(binomial_results$best_model)

# Red king crab

species_code <- 69322
common_name <- "red king crab"

# Setup data
dat <- cpue_1998[cpue_1998$SPECIES_CODE == species_code, ]
dat_bootstrap <- 
  draw_bootstrap_samples(
    x = dat, 
    seed = seed, 
    grouping_var = "TOW_PAIR", 
    n_draws = 1000, 
    replace = TRUE
  )
dat_oos <- cpue_other[cpue_other$SPECIES_CODE == species_code, ]

# Fit models
ols_results <- fit_ols(x = dat, bootstrap_samples = dat_bootstrap)
beta_results <- fit_betareg(x = dat, bootstrap_samples = dat_bootstrap)
binomial_results <- fit_binomial(x = dat, bootstrap_samples = dat_bootstrap)

ols_results$anderson_darling
ols_results$kurtosis
ols_results$cor_test
ols_results$bootstrap_results$ci
ols_results$fpc


p1 <- 
  check_ols_heteroskedasticity(
    dat = dat, 
    model = ols_results$best_model, 
    predictor = "CPUE_NO_KM2_15", 
    x_axis_name = expression(CPUE[15]*' (#/'*km^2*')')
  )

p2 <- 
  check_ols_heteroskedasticity(
    dat = dat, 
    model = ols_results$best_model, 
    predictor = "CPUE_NO_KM2_30", 
    x_axis_name = expression(CPUE[30]*' (#/'*km^2*')')
  )


cowplot::plot_grid(
  p1, p2
)



fpc_1998 <- 
  dplyr::bind_rows(
    binomial_results$bootstrap_results$ci,
    ols_results$bootstrap_results$ci,
    beta_results$bootstrap_results$ci
  )

fpc_table_1998 <- 
  fpc_1998 |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(fpc_table_1998)

loocv_table_1998 <- 
  dplyr::bind_rows(
    binomial_results$loocv_table,
    ols_results$loocv_table,
    beta_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(loocv_table_1998)

oos_table <- 
  indpendent_fpc_validation(
    fpc = fpc_1998,
    validation_dat = dat_oos
  ) |> 
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(oos_table)


rownames(fpc_table_1998) <- NULL
rownames(loocv_table_1998) <- NULL
rownames(oos_table) <- NULL

oos_table <- 
  oos_table |>
  dplyr::mutate(method = ifelse(is.na(method), model_name, method))


# Plot cross-year validation results

p_rmse_1998 <- 
  ggplot() +
  geom_point(
    data = oos_table,
    mapping = aes(x = method, y = rmse),
    size = rel(2.2)) +
  geom_errorbar(
    data = oos_table,
    mapping = aes(x = method, ymin = rmse_lci, ymax = rmse_uci),
    width = 0,
    linewidth = 1.05
  ) +
  scale_y_continuous(name = expression(RMSE*' (#/'*km^2*')')) +
  facet_wrap(~common_name, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )

p_tpe_1998 <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    data = oos_table,
    mapping = aes(x = common_name, y = tpe, color = method), position = position_dodge(width = 0.5),
    size = rel(2.2)) +
  geom_errorbar(
    data = oos_table,
    mapping = aes(x = common_name, ymin = tpe_lci, ymax = tpe_uci, color = method), position = position_dodge(width = 0.5),
    width = 0,
    linewidth = 1.05
  ) +
  scale_y_continuous(name = "Total percentage error (%)", limits = c(-100, 100), oob = oob_squish, expand = c(0, 0)) +
  scale_color_colorblind(name = "Method") +
  theme_bw() +
  theme(axis.title.x = element_blank())


loocv_table_1998_output <- 
  loocv_table_1998 |> 
  dplyr::mutate(
    years = 1998,
    sex_stage = "all"
  ) |>
  dplyr::select(common_name, years, model_name, method, rmse, tpe)

# Percentage difference (RMSE)
loocv_table_1998_output |> 
  dplyr::group_by(common_name, years) |>
  dplyr::summarise(min_rmse = min(rmse)) |>
  dplyr::inner_join(loocv_table_1998) |>
  dplyr::mutate(pdiff_rmse = (rmse-min_rmse)/min_rmse*100,
                pabs_diff_tpe = (abs(tpe)-min_rmse)/min_rmse*100)

fpc_table_1998_output <- 
  fpc_table_1998 |> 
  dplyr::mutate(
    years = 1998,
    sex_stage = "all"
  ) |>
  dplyr::mutate(
    ratio = paste0(round(ratio, 3), " (", round(ratio_lci, 3), "-", round(ratio_uci, 3), ")")
  ) |>
  dplyr::select(common_name, years, model_name, method, ratio)

xlsx::write.xlsx(
  loocv_table_1998_output, 
  file = here::here("analysis", "somerton_2002", "plots", "loocv_table_1998.xlsx"), 
  row.names = FALSE
)

xlsx::write.xlsx(
  fpc_table_1998_output, 
  file = here::here("analysis", "somerton_2002", "plots", "fpc_table_1998.xlsx"), 
  row.names = FALSE
)


# Predictions for the full data set

cpue_species <- 

fpc_table_1998$ratio[1] * cpue_other



# Fitting to the full data set ---------------------------------------------------------------------


