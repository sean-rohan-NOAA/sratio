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
    
    model_list <- 
      list(
        ols1 = ols1
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

# LOOCV models
run_loocv <- function(model_list, dat, mapper, ...) {
  
  results_list <- lapply(names(model_list), function(m_name) {
    mod <- model_list[[m_name]]
    n_obs <- nrow(dat)
    
    preds <- numeric(n_obs)
    
    for(jj in 1:n_obs) {
      # Update model excluding one observation
      fit_loocv <- update(mod, data = dat[-jj, , drop = FALSE])
      
      # Use mapper function to get the back-transformed prediction
      preds[jj] <- mapper(fit_loocv, dat[jj, , drop = FALSE], ...)
    }
    
    obs <- dat$CPUE_NO_KM2_30
    
    # Root mean square error
    rmse <- sqrt(mean((preds - obs)^2))
    
    # Total percentage error (relative bias in total count)
    tpe  <- 100 * (sum(preds * dat$AREA_SWEPT_KM2_30) - sum(dat$COUNT_30)) / sum(dat$COUNT_30)
    
    data.frame(model_name = m_name, rmse = rmse, tpe = tpe)
    
  })
  
  results <- do.call(rbind, results_list)
  results$best <- results$rmse == min(results$rmse)
  return(results)
}


# OLS ratio
ols_mapper <- function(model, test_row, bias_correct = TRUE) {
  pred_log_ratio <- predict(model, newdata = test_row)
  # Applying 0.5 * sigma2 bias correction
  sigma2 <- summary(model)$sigma^2
  
  if(bias_correct) { 
    ratio <- exp(pred_log_ratio + 0.5 * sigma2) 
    } else {
    ratio <- exp(pred_log_ratio)
  }
  
  return(test_row$CPUE_NO_KM2_15 / ratio)
}

# Binomial and beta
ccr_mapper <- function(model, test_row) {
  
  # Set combined count to 1 since it won't be observed but is required for predicting the response;
  # Although this doesn't have any effect on the predicted proportion since obs weights are not used
  # when type = "response"
  test_row$COMBINED_COUNT <- 1
  
  # Get predicted proportion (inv-logit scale)
  p <- predict(model, newdata = test_row, type = "response")
  ratio <- (1 - p) / p
  return(test_row$CPUE_NO_KM2_15 / ratio)
}

# Poisson and negative binomial
count_mapper <- function(model, test_row) {
  # Get predicted counts
  pred_count <- predict(model, newdata = test_row, type = "response")
  
  return(pred_count / test_row$AREA_SWEPT_KM2_30)
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
