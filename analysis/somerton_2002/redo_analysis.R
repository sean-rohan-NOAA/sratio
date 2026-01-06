# library(sratio)
library(nortest)
library(xlsx)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(dplyr)
library(here)

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
    YEAR = 1998,
    CPUE_NO_KM2_15 = COUNT_15 / AREA_SWEPT_KM2_15,
    CPUE_NO_KM2_30 = COUNT_30 / AREA_SWEPT_KM2_30,
    LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
    LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
    CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
    CPUE_RATIO = CPUE_NO_KM2_15/CPUE_NO_KM2_30,
    COMBINED_COUNT = ceiling(COUNT_30 + COUNT_15),
    PROP_15 = CPUE_NO_KM2_15/(CPUE_NO_KM2_15+CPUE_NO_KM2_30),
    EFFORT_RATIO = AREA_SWEPT_KM2_15/AREA_SWEPT_KM2_30,
    COUNT_15 = ceiling(COUNT_15),
    COUNT_30 = ceiling(COUNT_30),
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
  dplyr::summarise(COUNT = ceiling(sum(SAMPLING_FACTOR))) |>
  dplyr::inner_join(
    dplyr::select(
      sratio::data_1530$haul,
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
    COMBINED_COUNT = ceiling(COUNT_30 + COUNT_15),
    PROP_15 = CPUE_NO_KM2_15/(CPUE_NO_KM2_15+CPUE_NO_KM2_30),
    EFFORT_RATIO = AREA_SWEPT_KM2_15/AREA_SWEPT_KM2_30,
    common_name = sratio::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0)

saveRDS(object = cpue_other, file = here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))




# Load data ----

cpue_other <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))
cpue_1998 <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_1998.rds"))


# Model-fitting functions --------------------------------------------------------------------------

fit_ols_models <- 
  function(x) {

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
  
  aic_table <- 
    make_aic_table(
      model_list = model_list
  )
  
  loocv_table <- 
    rbind(
      run_loocv(model_list = model_list, dat = x, mapper = ols_mapper, bias_correct = TRUE),
      run_loocv(model_list = model_list, dat = x, mapper = ols_mapper, bias_correct = FALSE)
    )
  
  loocv_table$method <- c("OLS mean", "OLS median")
  
  fit_table <- 
    rbind(
      predict_fits(model_list = model_list, dat = x, mapper = ols_mapper, bias_correct = TRUE) |>
        dplyr::mutate(method = "OLS mean"),
      predict_fits(model_list = model_list, dat = x,  mapper = ols_mapper, bias_correct = FALSE) |>
        dplyr::mutate(method = "OLS median")
    )
  
  output <- list(
    models = model_list,
    best_model = model_list[[loocv_table$model_name[loocv_table$best][1]]],
    aic_table = aic_table,
    loocv_table = loocv_table,
    fit_table = fit_table
  )
  
  output[['cor_test']] <- list(
    CPUE_30 = summary(lm(abs(rstandard(output[['best_model']]))~I(log(x$CPUE_NO_KM2_30)))),
    CPUE_15 = summary(lm(abs(rstandard(output[['best_model']]))~I(log(x$CPUE_NO_KM2_15))))
  )
  
  output[['anderson_darling']] <- nortest::ad.test(rstandard(output[['best_model']]))
  output[['kurtosis']] <- e1071::kurtosis(rstandard(output[['best_model']]), type = 2)
  
  fpc <- miller_bias_correct(output[['best_model']])
  fpc$model_name <- loocv_table$model_name[loocv_table$best][1]
  
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

fit_lognormal <-
  function(x) {
    
    # Fit lognormal models
    lognormal_formulas <- 
      expand.grid(
        fixed_formula = c(
          CPUE_RATIO ~ 1,
          CPUE_RATIO ~ LOG_CPUE_NO_KM2_15,
          CPUE_RATIO ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
          CPUE_RATIO ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3),
          CPUE_RATIO ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3) + I(LOG_CPUE_NO_KM2_15^4)
        ),
        disp = NA
      )
    
    lognormal_index <- 1:nrow(lognormal_formulas)
    
    lognormal_models_list <- 
      lapply(lognormal_index, function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = lognormal_formulas$fixed_formula[[index]],
            family = lognormal(link = "identity"),
            weight = sqrt(COMBINED_COUNT),
            data = x
          )
        
        mod
        
      }
      
      )
    
    names(lognormal_models_list) <- paste0("lognormal", lognormal_index)

    # Only carry forward models that pass initial checks
    aic_table <- make_aic_table(lognormal_models_list)
    
    model_name_passed <- aic_table$model_name[aic_table$pass_check]
    
    model_list <- lognormal_models_list[names(lognormal_models_list) %in% model_name_passed]
    
    # Cross validation and fits
    loocv_table <- 
      run_loocv(model_list = model_list, dat = x, mapper = ratio_mapper)
    
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

fit_ccr_models <- 
  function(x) {
    
    # Fit beta catch comparison rate models
    ccr_beta_formulas <- 
      expand.grid(
        fixed_formula = c(
          PROP_15 ~ 1,
          PROP_15 ~ LOG_CPUE_NO_KM2_15,
          PROP_15 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
          PROP_15 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3),
          PROP_15 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3) + I(LOG_CPUE_NO_KM2_15^4)
        ),
        disp_formula =
          c(
            ~ 1,
            ~ LOG_CPUE_NO_KM2_15,
            ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2)
          )
      )
    
    ccr_beta_index <- 1:nrow(ccr_beta_formulas)
    
    ccr_beta_models_list <- 
      lapply(ccr_beta_index, function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = ccr_beta_formulas$fixed_formula[[index]],
            offset = log(EFFORT_RATIO),
            weight = sqrt(LOG_CPUE_NO_KM2_15),
            family = glmmTMB::beta_family(link = "logit"),
            disp = ccr_beta_formulas$disp_formula[[index]],
            data = x
          )
        
        mod
        
      }
      
      )
    
    names(ccr_beta_models_list) <- paste0("ccr_beta", ccr_beta_index)
    
    # Fit binomial models
    ccr_bin_formulas <- ccr_beta_formulas[, 1, drop = FALSE] |> unique()
    
    ccr_bin_index <- 1:nrow(ccr_bin_formulas)
    
    ccr_bin_models_list <- 
      lapply(ccr_bin_index, function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = ccr_bin_formulas$fixed_formula[[index]],
            offset = log(EFFORT_RATIO),
            weight = sqrt(LOG_CPUE_NO_KM2_15),
            family = binomial(link = "logit"),
            data = x
          )
        
        mod
        
      }
      
      )
    
    names(ccr_bin_models_list) <- paste0("ccr_bin", ccr_bin_index)
    
    
    # Only carry forward models that passed initial checks
    aic_table <- 
      dplyr::bind_rows(
        make_aic_table(
          ccr_beta_models_list
        ),
        make_aic_table(
          ccr_bin_models_list
        )
      )
    
    model_name_passed <- aic_table$model_name[aic_table$pass_check]
    
    model_list <- 
      c(
        ccr_beta_models_list,
        ccr_bin_models_list
      )
    
    model_list <- model_list[names(model_list) %in% model_name_passed]
    
    # Cross validation and fits
    loocv_table <-
      run_loocv(model_list = model_list, dat = x, mapper = ccr_mapper)
    
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

fit_prop_models <- 
  function(x) {
    
    # Fit beta-binomial models
    bb_formulas <- 
      expand.grid(
        fixed_formula = c(
          cbind(COUNT_15, COUNT_30) ~ 1,
          cbind(COUNT_15, COUNT_30) ~ LOG_CPUE_NO_KM2_15,
          cbind(COUNT_15, COUNT_30) ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
          cbind(COUNT_15, COUNT_30) ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3),
          cbind(COUNT_15, COUNT_30) ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3) + I(LOG_CPUE_NO_KM2_15^4)
        ),
        disp_formula =
          c(
            ~ 1,
            ~ LOG_CPUE_NO_KM2_15,
            ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2)
          )
      )
    
    bb_index <- 1:nrow(bb_formulas)
    
    bb_models_list <- 
      lapply(bb_index, function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = bb_formulas$fixed_formula[[index]],
            offset = log(EFFORT_RATIO),
            family = glmmTMB::betabinomial(link = "logit"),
            disp = bb_formulas$disp_formula[[index]],
            data = x
          )
        
        mod
        
      }
      
      )
    
    names(bb_models_list) <- paste0("bb", bb_index)
    
    # Fit binomial models
    bin_formulas <- bb_formulas[, 1, drop = FALSE] |> unique()
    
    bin_index <- 1:nrow(bin_formulas)
    
    bin_models_list <- 
      lapply(bin_index, function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = bin_formulas$fixed_formula[[index]],
            offset = log(EFFORT_RATIO),
            family = binomial(link = "logit"),
            data = x
          )
        
        mod
        
      }
      
      )
    
    names(bin_models_list) <- paste0("bin", bin_index)
    
    # Only carry forward models that passed initial checks
    aic_table <- 
      dplyr::bind_rows(
        make_aic_table(
          bb_models_list
        ),
        make_aic_table(
          bin_models_list
        )
      )
    
    model_name_passed <- aic_table$model_name[aic_table$pass_check]
    
    model_list <- 
      c(
        bb_models_list,
        bin_models_list
      )
    
    model_list <- model_list[names(model_list) %in% model_name_passed]
    
    loocv_table <-
      run_loocv(model_list = model_list, dat = x, mapper = prop_mapper)
    
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

fit_count_models <- 
  function(x) {
    
    # Fit negative binomial models
    nb_formulas <- 
      expand.grid(
        fixed_formula = c(
          COUNT_30 ~ 0 + LOG_CPUE_NO_KM2_15 ,
          COUNT_30 ~ 0 + LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
          COUNT_30 ~ 0 + LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3),
          COUNT_30 ~ 0 + LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3) + I(LOG_CPUE_NO_KM2_15^4),
          COUNT_30 ~ LOG_CPUE_NO_KM2_15,
          COUNT_30 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
          COUNT_30 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3),
          COUNT_30 ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2) + I(LOG_CPUE_NO_KM2_15^3) + I(LOG_CPUE_NO_KM2_15^4)
        ),
        disp_formula =
          c(
            ~ 1,
            ~ LOG_CPUE_NO_KM2_15,
            ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2)
          )
      )
    
    nb_index <- 1:nrow(nb_formulas)
    
    nb_models_list <- 
      lapply(nb_index, function(index) {
      
        mod <-
          glmmTMB::glmmTMB(
            formula = nb_formulas$fixed_formula[[index]],
            family = glmmTMB::nbinom1(link = "log"),
            offset = log(AREA_SWEPT_KM2_30),
            disp = nb_formulas$disp_formula[[index]],
            data = x
          )
        
        mod
        
      }
    
    )
    
    names(nb_models_list) <- paste0("nb", nb_index)
    
    # Fit Poisson models
    
    pois_formulas <- nb_formulas[, 1, drop = FALSE] |> unique()
    
    pois_index <- 1:nrow(pois_formulas)
    
    pois_models_list <- 
      lapply(pois_index, function(index) {
        
        mod <-
          glmmTMB::glmmTMB(
            formula = pois_formulas$fixed_formula[[index]],
            family = poisson(link = "log"),
            offset = log(AREA_SWEPT_KM2_30),
            data = x
          )
        
        mod
        
      }
      
      )
    
    names(pois_models_list) <- paste0("pois", pois_index)
    
    aic_table <- 
      dplyr::bind_rows(
        make_aic_table(
          pois_models_list
        ),
        make_aic_table(
          nb_models_list
        )
      )
    
    # Only carry forward models that passed initial checks
    model_name_passed <- aic_table$model_name[aic_table$pass_check]
    
    model_list <- 
      c(
        pois_models_list,
        nb_models_list
      )
    
    model_list <- model_list[names(model_list) %in% model_name_passed]
    
    # LOOCV and fits
    loocv_table <- 
      run_loocv(model_list = model_list, dat = x, mapper = count_mapper)
    
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
    pbias  <- 100 * (sum(preds - obs)) / sum(obs)

    data.frame(model_name = m_name, rmse = rmse, pbias = pbias)

  })
  
  results <- do.call(rbind, results_list)
  results$best <- results$rmse == min(results$rmse)
  return(results)
}

run_twofold_cv <- 
  function(model_list, validation_dat, mapper, subset = NULL, save_dir = NULL, common_name = NULL, ...) {
    
    obs_cpue    <- validation_dat$CPUE_NO_KM2_30
    obs_effort  <- validation_dat$AREA_SWEPT_KM2_30
    obs_count   <- validation_dat$COUNT_30
    
    results_list <- lapply(names(model_list), function(m_name) {
      mod <- model_list[[m_name]]
      
      preds <- mapper(
        model = mod,
        test_row = validation_dat,
        se_fit = TRUE,
        ...
      )
      
      # Root mean square error and total percentage error
      output <- 
        data.frame(
          model_name = m_name,
          rmse       = sqrt(mean((preds$fit - obs_cpue)^2)),
          rmse_lci   = sqrt(mean((preds$fit_lwr - obs_cpue)^2)),
          rmse_uci   = sqrt(mean((preds$fit_upr - obs_cpue)^2)),
          pbias        = 100 * sum(preds$fit - obs_cpue) / sum(obs_cpue),
          pbias_lci    = 100 * sum(preds$fit_lwr - obs_cpue) / sum(obs_cpue),
          pbias_uci    = 100 * sum(preds$fit_upr - obs_cpue) / sum(obs_cpue)
        )
      
      if(!is.null(save_dir)) {
        
        p1 <- ggplot() +
          geom_point(
            mapping = aes(x = obs_cpue, y = preds$fit)
          ) +
          geom_errorbar(
            mapping = aes(x = obs_cpue, ymin = preds$fit_lwr, ymax = preds$fit_upr),
            width = 0
          ) +
          geom_abline(slope = 1, intercept = 0, linetype = 2) +
          ggtitle(paste0(common_name, " ", m_name, " ", subset)) +
          scale_x_log10(name = expression("Observed CPUE[30]")) +
          scale_y_log10(name = expression("Predicted CPUE[30]"))
        
        fname <- paste0(
          "OBS_VS_PRED_", 
          gsub(x = common_name, pattern = " ", replacement = "_"), 
          "_", subset,
          "_", gsub(x = m_name, pattern = " ", replacement = "_"), 
          ".png")
        
        png(here::here(save_dir, fname), width = 140, height = 140, units = "mm", res = 300)
        print(p1)
        dev.off()
        
      }
      
      output
    })
    
    results <- do.call(what = rbind, args = results_list)
    results$best <- results$rmse == min(results$rmse)
    
    return(results)
    
  }

# Mapping functions for converting model predictions to estimation scale ---------------------------

# OLS ratio
ols_mapper <- function(model, test_row, bias_correct = TRUE, se_fit= FALSE) {
  
  sigma2 <- summary(model)$sigma^2
  
  if(se_fit) {
    pred_log_ratio <- predict(model, newdata = test_row, se.fit = TRUE)

    if(bias_correct) {
      ratio     <- exp(pred_log_ratio$fit + 0.5 * sigma2)
      ratio_lwr <- exp(pred_log_ratio$fit - 2 * pred_log_ratio$se.fit + 0.5 * sigma2)
      ratio_upr <- exp(pred_log_ratio$fit + 2 * pred_log_ratio$se.fit + 0.5 * sigma2)
    } else {
      ratio     <- exp(pred_log_ratio$fit)
      ratio_lwr <- exp(pred_log_ratio$fit - 2 * pred_log_ratio$se.fit)
      ratio_upr <- exp(pred_log_ratio$fit + 2 * pred_log_ratio$se.fit)
    }
    
    output <- 
      data.frame(
        fit     = test_row$CPUE_NO_KM2_15 / ratio,
        fit_lwr = test_row$CPUE_NO_KM2_15 / ratio_lwr,
        fit_upr = test_row$CPUE_NO_KM2_15 / ratio_upr
      )
    
  } else {
    pred_log_ratio <- predict(model, newdata = test_row)
    
    if(bias_correct) { 
      # Applying 0.5 * sigma2 bias correction
      ratio <- exp(pred_log_ratio + 0.5 * sigma2) 
    } else {
      ratio <- exp(pred_log_ratio)
    }
    
    output <- test_row$CPUE_NO_KM2_15 / ratio
  }
  
  return(output)
  
}

# lognormal ratio
ratio_mapper <- function(model, test_row, se_fit = FALSE) {
  # Get predicted ratio
  
  if(se_fit) {
    # Two-fold CV
    inv_link  <- family(model)$linkinv
    
    ratio     <- predict(model, newdata = test_row, type = "link", se.fit = TRUE)
    ratio_lwr <- ratio$fit - 2* ratio$se.fit
    ratio_upr <- ratio$fit + 2* ratio$se.fit 
    
    ratio     <- inv_link(ratio$fit)
    ratio_lwr <- inv_link(ratio_lwr)
    ratio_upr <- inv_link(ratio_upr)
    
    output <- 
      data.frame(
      fit =     test_row$CPUE_NO_KM2_15 / ratio,
      fit_lwr = test_row$CPUE_NO_KM2_15 / ratio_lwr,
      fit_upr = test_row$CPUE_NO_KM2_15 / ratio_upr
    )
    
  } else {
    # LOOCV
    ratio <- predict(model, newdata = test_row, type = "response")
    output <- test_row$CPUE_NO_KM2_15 / ratio
  }
  
  return(output)
}

# Binomial and beta-binomial catch comparison rate models
ccr_mapper <- function(model, test_row, se_fit = FALSE) {

  if(se_fit) {
    # Two-fold CV
    inv_link  <- family(model)$linkinv

    fit         <- predict(model, newdata = test_row, type = "link", se.fit = TRUE)
    fit_lwr     <- fit$fit - 2* fit$se.fit
    fit_upr     <- fit$fit + 2* fit$se.fit
    
    # Response scale (CCR)
    ccr     <- inv_link(fit$fit)
    ccr_lwr <- inv_link(fit_lwr)
    ccr_upr <- inv_link(fit_upr)

    # Fishing power
    ratio     <-  ccr / (1 - ccr)
    ratio_lwr <-  ccr_lwr / (1 - ccr_lwr)
    ratio_upr <-  ccr_upr / (1 - ccr_upr)



    output <-
      data.frame(
        fit =     test_row$CPUE_NO_KM2_15 / ratio,
        fit_lwr = test_row$CPUE_NO_KM2_15 / ratio_lwr,
        fit_upr = test_row$CPUE_NO_KM2_15 / ratio_upr
      )

  } else {
    # LOOCV; Get predicted proportion (inv-logit scale)
    ccr <- predict(model, newdata = test_row, type = "response")
    ratio <- ccr / (1 - ccr)
    output <- test_row$CPUE_NO_KM2_15 / ratio
  }

  return(output)
}


prop_mapper <- function(model, test_row, se_fit = FALSE) {
  
  if(se_fit) {
    # Two-fold CV
    inv_link  <- family(model)$linkinv
    
    rr         <- predict(model, newdata = test_row, type = "link", se.fit = TRUE)
    rr_lwr     <- rr$fit - 2* rr$se.fit
    rr_upr     <- rr$fit + 2* rr$se.fit 
    
    p     <- inv_link(rr$fit)
    p_lwr <- inv_link(rr_lwr)
    p_upr <- inv_link(rr_upr)
    
    output <- 
      data.frame(
        fit =     (test_row$COUNT_15 / p - test_row$COUNT_15) / test_row$AREA_SWEPT_KM2_30,
        fit_lwr = (test_row$COUNT_15 / p_lwr - test_row$COUNT_15) / test_row$AREA_SWEPT_KM2_30,
        fit_upr = (test_row$COUNT_15 / p_upr - test_row$COUNT_15) / test_row$AREA_SWEPT_KM2_30
      )
    
  } else {
    # LOOCV; Get predicted proportion (inv-logit scale)
    p <- predict(model, newdata = test_row, type = "response")
    output <- (test_row$COUNT_15 / p - test_row$COUNT_15) / test_row$AREA_SWEPT_KM2_30
  }
  
  return(output)
  
}

# Poisson and negative binomial
count_mapper <- function(model, test_row, se_fit = FALSE) {
  
  if(se_fit) {
    
    inv_link  <- family(model)$linkinv
    
    mu     <- predict(model, newdata = test_row, type = "link", se.fit = TRUE)
    mu_lwr <- mu$fit - 2 * mu$se.fit
    mu_upr <- mu$fit + 2 * mu$se.fit 
    
    pred_count     <- inv_link(mu$fit)
    pred_count_lwr <- inv_link(mu_lwr)
    pred_count_upr <- inv_link(mu_upr)
    
    output <- 
      data.frame(
        fit =     pred_count / test_row$AREA_SWEPT_KM2_30,
        fit_lwr = pred_count_lwr / test_row$AREA_SWEPT_KM2_30,
        fit_upr = pred_count_upr / test_row$AREA_SWEPT_KM2_30
      )
    
  } else {
    # LOOCV; Get predicted counts
    pred_count <- predict(model, newdata = test_row, type = "response")
    output <- pred_count / test_row$AREA_SWEPT_KM2_30
  }

  
  return(output)
}

# Make predictions ---------------------------------------------------------------------------------

predict_fits <- 
  function(model_list, dat, mapper, ...) {
    
    log_cpue <- seq(log(min(dat$CPUE_NO_KM2_15)), log(max(dat$CPUE_NO_KM2_15)), length.out = 300)
    
    fit_data <-
      data.frame(
        AREA_SWEPT_KM2_15 = 0.02,
        AREA_SWEPT_KM2_30 = 0.04,
        EFFORT_RATIO = 0.02/0.04,
        LOG_CPUE_NO_KM2_15 = log_cpue,
        CPUE_NO_KM2_15  = exp(log_cpue),
        COMBINED_COUNT = 1
      )
    
    fit_data$COUNT_15 <- fit_data$CPUE_NO_KM2_15*fit_data$AREA_SWEPT_KM2_15
    
    results_list <- 
      lapply(
        names(model_list), 
        function(m_name) {
          
          output <- 
            cbind(
              fit_data,
              # mapper(model = model_list[[m_name]], test_row = fit_data, se_fit = TRUE) testing
              mapper(model = model_list[[m_name]], test_row = fit_data, se_fit = TRUE, ...)
            )
          
          output$model_name <- m_name
          
          output
          
        }
        
      )
    
    return(do.call(what = rbind, args = results_list))
    
  }



# Diagnostics --------------------------------------------------------------------------------------
make_aic_table <- 
  function(model_list) {
    
    results <- data.frame(
      model_name = names(model_list) %||% seq_along(model_list),
      formula    = sapply(model_list, function(m) paste(format(formula(m)), collapse = "")),
      
      # Dispersion formula from glmmTMB
      disp       = sapply(
        model_list, 
        function(m){
          if(is(m, "lm")) {
            out <- NA} else{
              out <- paste(format(m$modelInfo$allForm$dispformula))
            }
          out
        }),
      aic        = round(sapply(model_list, AIC), 2),
      k          = sapply(model_list, function(m) attr(logLik(m), "df")),
      convergence = sapply(
        model_list, 
        function(m) {
          if(is(m, "lm")) {
            out <- NA} else{
              out <- m$fit$convergence
            }
          out
        }),
      pdhess = sapply(
        model_list, 
        function(m) {
          if(is(m, "lm")) {
            out <- NA} else{
              out <- m$sdr$pdHess
            }
          out
        }),
      max_gradient = sapply(
        model_list, 
        function(m) {
          if(is(m, "lm")) {
            out <- NA} else{
              out <- max(abs(m$sdr$gradient.fixed))
            }
          out
        }),
      stringsAsFactors = FALSE
    )
    
    results$pass_check <- 
      results$convergence == 0 & results$pdhess & abs(results$max_gradient) < 0.001
    
    results$pass_check[grepl(pattern = "ols", x = results$model_name)] <- TRUE
    
    results$delta_aic <- results$aic - min(results$aic, na.rm = TRUE)
    
    candidates <- results[results$delta_aic < 2, ]
    
    return(results[order(results$aic), ])
    
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

dharma_plots <- 
  function(model_list, common_name, subset = NULL, nsim = 1000, save_dir) {
    
    model_names    <- names(model_list)
    results        <- vector(mode = "list", length(model_list))
    names(results) <- model_names
    
    for(ii in seq_along(model_list)) {
      
      results[[ii]] <- try(DHARMa::simulateResiduals(model_list[[ii]], n = 1000), silent = TRUE)
      
      if(is(results[[ii]], "try-error")) {
        warning(paste0("warning: Unable to produce DHARMa residuals for model ", model_names[ii]))
        next
      }
      
      fname <- paste0(
        "DHARMa_", 
        gsub(x = common_name, pattern = " ", replacement = "_"), 
        "_", subset,
        "_", gsub(x = model_names[ii], pattern = " ", replacement = "_"), 
        ".png")
      
      png(here::here(save_dir, fname), width = 169, height = 120, units = "mm", res = 300)
      print(DHARMa::plotResiduals(results[[ii]]))
      dev.off()
      
    }
    
    return(results)
    
  }

# Bootstrap function -------------------------------------------------------------------------------
draw_bootstrap_samples <- 
  function(x, seed = NULL, grouping_var, n_draws = 1000, replace = TRUE) {
    
    grouping_var_values <- unique(x[[grouping_var]])
    n_obs <- length(grouping_var_values)
    
    draws <- vector(mode = "list", length = n_draws)
    
    if(is.numeric(seed)) {
      set.seed(seed)
    }
    
    for(ii in 1:n_draws) {
      
      samp <- 
        data.frame(
          var = sample(x = grouping_var_values, replace = replace, size = n_obs)
          )
      
      draws[[ii]] <-
        dplyr::inner_join(
          x, samp, 
          by = setNames("var", grouping_var), 
          relationship = "many-to-many"
      ) |>
        dplyr::mutate(draw = ii)
    }
    
    return(draws)
    
  }

# Wrapper function to run analyses -----------------------------------------------------------------
run_analysis <- 
  function(dat, dat_oos = NULL, common_name, subset_name, contrast_name = NULL) {
    
    fits_dir <- here::here("analysis", "somerton_2002", "plots", paste0(subset_name, "_fits"))
    dharma_dir <- here::here("analysis", "somerton_2002", "plots", paste0(subset_name, "_fits"), "dharma")
    obs_pred_dir <- here::here("analysis", "somerton_2002", "plots", paste0(subset_name, "_fits"), "obs_pred_dir")
    dir.create(fits_dir, showWarnings = FALSE)
    dir.create(dharma_dir, showWarnings = FALSE)
    dir.create(obs_pred_dir, showWarnings = FALSE)
    
    
    # Fit models, check convergence
    ols_results <- fit_ols_models(x = dat)
    ccr_results <- fit_ccr_models(x = dat)
    prop_results <- fit_prop_models(x = dat)
    lognormal_results <- fit_lognormal(x = dat)
    count_results <- fit_count_models(x = dat)
    
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
          ) |>
            dplyr::mutate(method = "OLS mean"),
          run_twofold_cv(
            model_list = ols_results$models, 
            validation_dat = dat_oos, 
            mapper = ols_mapper, 
            common_name = common_name,
            subset = contrast_name,
            save_dir = obs_pred_dir,
            bias_correct = FALSE
          ) |>
            dplyr::mutate(method = "OLS median"),
          run_twofold_cv(
            model_list = prop_results$models, 
            validation_dat = dat_oos, 
            mapper = prop_mapper, 
            common_name = common_name,
            subset = contrast_name,
            save_dir = obs_pred_dir
          )
          )|>
        dplyr::mutate(common_name = common_name) |> 
        dplyr::arrange(rmse) |>
        dplyr::mutate(method = ifelse(is.na(method), model_name, method))
      
    }
    
    # Check DHARMa residuals
    dharma_results <- 
      dharma_plots(
        model_list = c(
          ols_results$models,
          ccr_results$models,
          lognormal_results$models,
          prop_results$models,
          count_results$models
        ),
        common_name = common_name,
        nsim = 1000,
        subset = subset_name,
        save_dir = dharma_dir
      )
    
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
    
    p_heteroskedasticity <- cowplot::plot_grid(p1, p2)
    
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
        dharma_results = dharma_results,
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

# Fit 2021 to 2024 (2CV: 1998 + 1995 2CV) ----------------------------------------------------------
# results_2124 <- vector(mode = "list", length = nrow(set_species_all))
# names(results_2124) <- set_species_all$common_name
# 
# for(ll in 1:nrow(set_species_all)) {
# 
#   species_code <- set_species_all$species_code[ll]
#   sex <- set_species_all$sex[ll]
# 
#   sel_dat <-
#     cpue_other[cpue_other$SPECIES_CODE == species_code & cpue_other$SEX == sex & cpue_other$YEAR != 1995, ]
# 
#   sel_oos <-
#     dplyr::bind_rows(
#       cpue_other[cpue_other$SPECIES_CODE == species_code & cpue_other$SEX == sex & cpue_other$YEAR == 1995, ],
#       cpue_1998[cpue_1998$SPECIES_CODE == species_code & cpue_1998$SEX == sex, ]
#     )
# 
# 
#   results_2124[[ll]] <-
#     run_analysis(
#       dat = sel_dat,
#       dat_oos = sel_oos,
#       common_name = set_species_all$common_name[ll],
#       subset_name = "21to24",
#       contrast_name = "21to24FitVs95to98"
#     )
# 
# }
# 
# save(results_2124, file = here::here("analysis", "somerton_2002", "output", "results_21to24.rda"))



# Plot cross-year validation results ------------

p_ols <-
  ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    data = oos_table_ols,
             mapping = aes(x = common_name, y = pbias, color = method), position = position_dodge(width = 0.5)
    ) +
    geom_errorbar(
      data = oos_table_ols,
      mapping = aes(
        x = common_name, ymin = pbias_lci, ymax = pbias_uci, color = method
      ),
      position = position_dodge(width = 0.5),
      width = 0
    ) +
    scale_y_continuous(name = "Total percentage error (%)") +
    scale_color_colorblind(name = "Method") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.22, 0.82),
          axis.title.y = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))

png(
  filename =  here::here(fits_dir, "PBIAS_OLS_other_years_from_1998.png"),
  width = 80,
  height = 80, 
  res = 300, 
  units = "mm"
)
print(p_ols_1998)
dev.off()

p_rmse_1998 <-
  ggplot() +
  geom_point(
    data = oos_table_ols_1998,
    mapping = aes(x = method, y = rmse),
    size = rel(2.2)) +
  geom_errorbar(
    data = oos_table_ols_1998,
    mapping = aes(x = method, ymin = rmse_lci, ymax = rmse_uci),
    width = 0,
    linewidth = 1.05
  ) +
  scale_y_continuous(name = expression(RMSE*' (#/'*km^2*')')) +
  facet_wrap(~common_name, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9)
  )

png(
  filename = here::here("analysis", "somerton_2002", "plots", "1998_fits", "RMSE_OLS_other_years_from_1998.png"),
  width = 169,
  height = 60, 
  res = 300, 
  units = "mm"
)
print(p_rmse_1998)
dev.off()


test <- results_all[[4]]$fit_table


ggplot() +
  geom_path(
    data = dplyr::mutate(
      test, method = ifelse(is.na(method), model_name, method)) |>
      dplyr::inner_join(response_type),
    mapping = aes(x = CPUE_NO_KM2_15, y = fit, color = type)
  ) +
  geom_abline(slope =1 , intercept = 0, linetype = 2) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_colorblind() +
  facet_wrap(~method)

