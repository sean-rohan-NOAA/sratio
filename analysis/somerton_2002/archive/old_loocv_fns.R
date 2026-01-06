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