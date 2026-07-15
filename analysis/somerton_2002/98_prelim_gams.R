# Didn't use GAMs because they led to unstable predictions for large values

# GAM functions --------------------------------------------------------------------------

fit_gam_lognormal <-
  function(x) {
    
    lognormal1 <- 
      gamlss::gamlss(
        formula = CPUE_RATIO ~ 1,
        family = gamlss.dist::LOGNO2(),
        data = x
      )
    
    lognormal2 <- 
      gamlss::gamlss(
        formula = CPUE_RATIO ~ pb(LOG_CPUE_NO_KM2_15),
        family = gamlss.dist::LOGNO2(),
        data = x
      )
      
    model_list <- 
      list(
        lognormal1 = lognormal1,
        lognormal2 = lognormal2
      )
    
    aic_table <- 
      make_aic_table(
        model_list = model_list
      )
    
    loocv_table <- 
      run_gam_loocv(model_list = model_list, dat = x, mapper = ratio_mapper)
    
    fit_table <- 
      predict_fits(model_list = model_list, dat = x, mapper = ratio_mapper)
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      fit_table = fit_table
    )
    
    return(output)
    
  }

fit_gam_ccr_models <- 
  function(x) {
    
    ccr_beta1 <-
        gamlss::gamlss(
        formula = PROP_15 ~ 1,
        family = gamlss.dist::BE(mu.link = "logit", sigma.link = "logit"),
        weight = sqrt(CPUE_NO_KM2_15),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    ccr_beta2 <-
      gamlss::gamlss(
        formula = PROP_15 ~ 1,
        sigma.formula =  ~ pb(LOG_CPUE_NO_KM2_15),
        family = gamlss.dist::BE(mu.link = "logit", sigma.link = "logit"),
        weight = sqrt(CPUE_NO_KM2_15),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    ccr_beta3 <-
      gamlss::gamlss(
        formula = PROP_15 ~ pb(LOG_CPUE_NO_KM2_15),
        family = gamlss.dist::BE(mu.link = "logit", sigma.link = "logit"),
        weight = sqrt(CPUE_NO_KM2_15),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    ccr_beta4 <-
      gamlss::gamlss(
        formula = PROP_15 ~ pb(LOG_CPUE_NO_KM2_15),
        sigma.formula =  ~ pb(LOG_CPUE_NO_KM2_15),
        family = gamlss.dist::BE(mu.link = "logit", sigma.link = "logit"),
        weight = sqrt(CPUE_NO_KM2_15),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    ccr_bin1 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COMBINED_COUNT) ~ 1,
        family = gamlss.dist::BI(mu.link = "logit"),
        control = gamlss.control(trace = FALSE),
        data = x
      )

    
    ccr_bin2 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COMBINED_COUNT) ~ pb(LOG_CPUE_NO_KM2_15),
        family = gamlss.dist::BI(mu.link = "logit"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    model_list <- 
      list(
        ccr_beta1 = ccr_beta1,
        ccr_beta2 = ccr_beta2,
        ccr_beta3 = ccr_beta3,
        ccr_beta4 = ccr_beta4,
        ccr_bin1 = ccr_bin1,
        ccr_bin2 = ccr_bin2
      )
    
    aic_table <- 
      dplyr::bind_rows(
        make_aic_table(
          model_list[grepl(x = names(model_list), pattern = "beta")] 
        ),
        make_aic_table(
          model_list[grepl(x = names(model_list), pattern = "bin")] 
        )
      )
    
    loocv_table <-
      run_gam_loocv(model_list = model_list, dat = x, mapper = ccr_mapper)
    
    fit_table <- 
      predict_fits(model_list = model_list, dat = x, mapper = ccr_mapper)
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      fit_table = fit_table
    )
    
    # Residual plots
    
    return(output)
    
  }

fit_gam_prop_models <- 
  function(x) {
    
    bb1 <-
       gamlss::gamlss(
        formula = cbind(COUNT_15, COUNT_30) ~ 1,
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::BB(mu.link = "logit", sigma.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    bb2 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COUNT_30) ~ 1,
        sigma.formula = ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::BB(mu.link = "logit", sigma.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    bb3 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COUNT_30) ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::BB(mu.link = "logit", sigma.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    bb4 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COUNT_30) ~ pb(LOG_CPUE_NO_KM2_15),
        sigma.formula = ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::BB(mu.link = "logit", sigma.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    
    bin1 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COUNT_30) ~ 1,
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::BB(mu.link = "logit"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    bin2 <-
      gamlss::gamlss(
        formula = cbind(COUNT_15, COUNT_30) ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::BB(mu.link = "logit"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    
    model_list <- 
      list(
        bb1 = bb1,
        bb2 = bb2,
        bb3 = bb3,
        bb4 = bb4,
        bin1 = bin1,
        bin2 = bin2
      )
    
    aic_table <- 
      make_aic_table(
        model_list 
      )
    
    loocv_table <-
      run_gam_loocv(model_list = model_list, dat = x, mapper = prop_mapper)
    
    fit_table <- 
      predict_fits(model_list = model_list, dat = x, mapper = prop_mapper)
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      fit_table = fit_table
    )
    
    # Residual plots
    
    return(output)
    
  }

fit_gam_count_models <- 
  function(x) {
    
    pois1 <-
      gamlss::gamlss(
        formula = COUNT_30 ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::PO(mu.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
    
    nb1 <-
      gamlss::gamlss(
        formula = COUNT_30 ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::NBI(mu.link = "log", sigma.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
      
      nb2 <- gamlss::gamlss(
        formula = COUNT_30 ~ pb(LOG_CPUE_NO_KM2_15),
        sigma.formula = COUNT_30 ~ pb(LOG_CPUE_NO_KM2_15),
        offset = log(EFFORT_RATIO),
        family = gamlss.dist::NBI(mu.link = "log", sigma.link = "log"),
        control = gamlss.control(trace = FALSE),
        data = x
      )
      
    
    model_list <- 
      list(
        pois1 = pois1,
        nb2 = nb2
      )
    
    aic_table <- 
      dplyr::bind_rows(
        make_aic_table(
          model_list[grepl(pattern = "pois", x = names(model_list))]
        ),
        make_aic_table(
          model_list[grepl(pattern = "nb", x = names(model_list))]
        )
      )
    
    loocv_table <- 
      run_gam_loocv(model_list = model_list, dat = x, mapper = count_mapper)
    
    fit_table <- 
      predict_fits(model_list = model_list, dat = x, mapper = count_mapper)
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      fit_table = fit_table
    )
    
    return(output)
    
  }

# Cross validation functions -----------------------------------------------------------------------
run_gam_loocv <- function(model_list, dat, mapper, ...) {
  
  n_mod <- length(model_list)
  n_obs <- nrow(dat)
  results_list <- vector(mode = "list", length = n_mod)
  
  for(ii in 1:n_mod) {
    
    m_name <- names(model_list)[[ii]]
    
    mod <- model_list[[m_name]]
    
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
    pbias  <- 100 * (sum(preds - obs)) / sum(obs)
    
    results_list[[ii]] <- data.frame(model_name = m_name, rmse = rmse, pbias = pbias)
    
  }
  
  results <- do.call(rbind, results_list)
  results$best <- results$rmse == min(results$rmse)
  return(results)
}

# Wrapper function to run analyses -----------------------------------------------------------------
run_gam_analysis <- 
  function(dat, dat_oos = NULL, common_name, subset_name, contrast_name = NULL) {
    
    fits_dir <- here::here("analysis", "somerton_2002", "plots", paste0(subset_name, "_gamlss_fits"))
    obs_pred_dir <- here::here("analysis", "somerton_2002", "plots", paste0(subset_name, "_gamlss_fits"), "obs_pred_dir")
    dir.create(fits_dir, showWarnings = FALSE)
    dir.create(obs_pred_dir, showWarnings = FALSE)
    
    
    # Fit models, check convergence
    ols_results <- fit_ols_models(x = dat)
    ccr_results <- fit_gam_ccr_models(x = dat)
    prop_results <- fit_gam_prop_models(x = dat)
    lognormal_results <- fit_gam_lognormal(x = dat)
    count_results <- fit_gam_count_models(x = dat)
    
    # Two-fold CV 
    
    oos_table <- NA
    
    if(!is.null(dat_oos)) {
      
      oos_table <- 
        dplyr::bind_rows(
          run_twofold_cv(
            model_list = ccr_results$models, 
            validation_dat = dat_oos, 
            mapper = ccr_mapper, 
            common_name = common_name,
            subset = contrast_name,
            save_dir = obs_pred_dir
          ),
          run_twofold_cv(
            model_list = lognormal_results$models, 
            validation_dat = dat_oos, 
            mapper = ratio_mapper, 
            common_name = common_name,
            subset = contrast_name,
            save_dir = obs_pred_dir
          ),
          run_twofold_cv(
            model_list = count_results$models, 
            validation_dat = dat_oos, 
            mapper = count_mapper, 
            common_name = common_name,
            subset = contrast_name,
            save_dir = obs_pred_dir
          ),
          run_twofold_cv(
            model_list = ols_results$models, 
            validation_dat = dat_oos, 
            mapper = ols_mapper, 
            common_name = common_name,
            subset = contrast_name,
            save_dir = obs_pred_dir,
            bias_correct = TRUE
          )
        ) |>
        dplyr::mutate(common_name = common_name) |> 
        dplyr::arrange(rmse) |>
        dplyr::mutate(method = ifelse(is.na(method), model_name, method))
      
    }
    
    # Make AIC table that includes convergence checks
    aic_table <-
      dplyr::bind_rows(
        ols_results$aic_table,
        ccr_results$aic_table,
        prop_results$aic_table,
        lognormal_results$aic_table,
        count_results$aic_table
      ) |>
      dplyr::mutate(common_name = common_name)
    
    # Append cross validation results
    loocv_table <- 
      dplyr::bind_rows(
        ols_results$loocv_table,
        ccr_results$loocv_table,
        prop_results$loocv_table,
        lognormal_results$loocv_table,
        count_results$loocv_table
      ) |>
      dplyr::arrange(rmse) |>
      dplyr::mutate(common_name = common_name) |> 
      dplyr::mutate(method = ifelse(is.na(method), model_name, method))
    
    fit_table <-
      dplyr::bind_rows(
        ols_results$fit_table,
        ccr_results$fit_table,
        prop_results$fit_table,
        lognormal_results$fit_table,
        count_results$fit_table
      ) |>
      dplyr::mutate(
        common_name = common_name,
        subset_name = subset_name
      )
    
    converged_table <-
      aic_table |>
      dplyr::select(model_name, common_name, pass_check) |>
      dplyr::inner_join(
        loocv_table
      )
    
    output <- 
      list(
        dat = dat,
        dat_oos = dat_oos,
        contrast_name = contrast_name,
        subset_name = subset_name,
        common_name = common_name,
        ols_results = ols_results, 
        ccr_results = ccr_results, 
        prop_results = prop_results,
        lognormal_results = lognormal_results, 
        count_results = count_results,
        loocv_table = loocv_table,
        aic_table = aic_table,
        oos_table = oos_table,
        fit_table = fit_table,
        converged_table = converged_table,
        p_heteroskedasticity = p_heteroskedasticity
      )
    
    return(output)
    
  }

# Fitting Somerton et al. (2002) 1998 data ---------------------------------------------------------

set_species_2002 <- 
  data.frame(
    species_code = c(68560, 68580, 69322),
    common_name = c("Tanner crab", "snow crab", "red king crab")
  )

results_2002 <- vector(mode = "list", length = nrow(set_species_2002))
names(results_2002) <- set_species_2002$common_name

for(kk in 1:nrow(set_species_2002)) {
  
  species_code <- set_species_2002$species_code[kk]
  
  results_2002[[kk]] <- 
    run_analysis(
      dat = cpue_1998[cpue_1998$SPECIES_CODE == species_code, ], 
      dat_oos = cpue_other[cpue_other$SPECIES_CODE == species_code, ], 
      common_name = set_species_2002$common_name[kk], 
      subset_name = "2002", 
      contrast_name = "2002FitVsOther"
    )
  
}

save(results_2002, file = here::here("analysis", "somerton_2002", "output", "results_2002.rda"))


# All years (no twofold CV) ------------------------------------------------------------------------
set_species_all <- 
  data.frame(
    species_code = c(68560, 68560, 68580, 68580, 69322, 69322),
    sex = c("M", "F", "M", "F", "M", "F"),
    common_name = 
      c(
        "Tanner crab (male)", 
        "Tanner crab (female)", 
        "snow crab (male)", 
        "snow crab (female)", 
        "red king crab (male)", 
        "red king crab (female)"
      )
  )

results_all <- vector(mode = "list", length = nrow(set_species_all))
names(results_all) <- set_species_all$common_name

for(kk in 1:nrow(set_species_all)) {
  
  species_code <- set_species_all$species_code[kk]
  sex <- set_species_all$sex[kk]
  
  sel_dat <- dplyr::bind_rows(
    cpue_1998[cpue_1998$SPECIES_CODE == species_code & cpue_1998$SEX == sex, ],
    cpue_other[cpue_other$SPECIES_CODE == species_code & cpue_other$SEX == sex, ]
  )
  
  results_all[[kk]] <- 
    run_analysis(
      dat = sel_dat, 
      common_name = set_species_all$common_name[kk], 
      subset_name = "All"
    )
  
}

save(results_all, file = here::here("analysis", "somerton_2002", "output", "results_all.rda"))


# Fit 1998 (2CV: other years 2CV) ------------------------------------------------------------------
results_1998 <- vector(mode = "list", length = nrow(set_species_all))
names(results_1998) <- set_species_all$common_name

for(ll in 1:nrow(set_species_all)) {
  
  species_code <- set_species_all$species_code[ll]
  sex <- set_species_all$sex[ll]
  
  results_1998[[ll]] <- 
    run_analysis(
      dat = cpue_1998[cpue_1998$SPECIES_CODE == species_code & cpue_1998$SEX == sex, ], 
      dat_oos = cpue_other[cpue_other$SPECIES_CODE == species_code & cpue_other$SEX == sex, ], 
      common_name = set_species_all$common_name[ll], 
      subset_name = 1998,
      contrast_name = "1998FitVsOther"
    )
  
}

save(results_1998, file = here::here("analysis", "somerton_2002", "output", "results_1998.rda"))


# Fit 1995 + 2021 to 2024 (2CV: 1998 2CV) ----------------------------------------------------------
results_95_2124 <- vector(mode = "list", length = nrow(set_species_all))
names(results_95_2124) <- set_species_all$common_name

for(ll in 1:nrow(set_species_all)) {
  
  species_code <- set_species_all$species_code[ll]
  sex <- set_species_all$sex[ll]
  
  results_95_2124[[ll]] <- 
    run_analysis(
      dat = cpue_other[cpue_other$SPECIES_CODE == species_code & cpue_other$SEX == sex, ], 
      dat_oos = cpue_1998[cpue_1998$SPECIES_CODE == species_code & cpue_1998$SEX == sex, ], 
      common_name = set_species_all$common_name[ll], 
      subset_name = "95_2124",
      contrast_name = "95_2124FitVs1998"
    )
  
}

save(results_95_2124, file = here::here("analysis", "somerton_2002", "output", "results_95_2124.rda"))




test <- expand.grid(
  fixed_formula = c(
    COUNT_30 ~ 0 + LOG_CPUE_NO_KM2_15 ,
    COUNT_30 ~ 0 + LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
    COUNT_30 ~ LOG_CPUE_NO_KM2_15,
    COUNT_30 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2)
  ),
  sigma_formula =
    c(
      ~ 1,
      ~ LOG_CPUE_NO_KM2_15,
      ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2)
    )
)

test[1,]
test$fixed_formula[1]



