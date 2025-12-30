library(sratio)
library(nortest)
library(xlsx)
library(glmmTMB)
library(DHARMa)



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
  

dharma_dir <- here::here("analysis", "somerton_2002", "plots", "dharma")
dir.create(dharma_dir, showWarnings = FALSE)

# Load data ----

cpue_other <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))
cpue_1998 <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_1998.rds"))


# CPUE ratio regression ----- 

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
              out <- paste(format(m$modlInfo$allForm$dispformula))
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
    
    results$pass_check <- results$convergence == 0 & results$pdhess & abs(results$max_gradient) < 0.001
    
    results$delta_aic <- results$aic - min(results$aic, na.rm = TRUE)
    
    candidates <- results[results$delta_aic < 2, ]
    
    return(results[order(results$aic), ])
  
} 



fit_ols_models <- 
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
    rbind(
      run_loocv(model_list = model_list, dat = dat, mapper = ols_mapper, bias_correct = TRUE),
      run_loocv(model_list = model_list, dat = dat, mapper = ols_mapper, bias_correct = FALSE)
    )
  
  loocv_table$method <- c("OLS mean", "OLS median")
  
  output <- list(
    models = model_list,
    best_model = model_list[[loocv_table$model_name[loocv_table$best][1]]],
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


fit_lognormal <-
  function(x, bootstrap_samples = NULL) {
    
    lognormal1 <- 
      glmmTMB::glmmTMB(
        formula = CPUE_RATIO ~ 1,
        family = lognormal(link = "log"),
        weight = sqrt(COMBINED_COUNT),
        data = x
      )
    
    model_list <- 
      list(
        lognormal1 = lognormal1
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
        run_loocv(model_list = model_list, dat = dat, mapper = ratio_mapper)
    
    output <- list(
      models = model_list,
      best_model = model_list[[loocv_table$model_name[loocv_table$best]]],
      aic_table = aic_table,
      loocv_table = loocv_table,
      bootstrap_results = bootstrap_results
    )
    
    return(output)
    
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



fit_ccr_models <- 
  function(x, bootstrap_samples = NULL) {
    
    bb1 <- 
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        family = betabinomial(link = "logit"),
        offset = log(EFFORT_RATIO),
        weight = COMBINED_COUNT,
        data = x
      )
    
    bb2 <- 
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        dispformula = ~ LOG_CPUE_NO_KM2_15,
        family = betabinomial(link = "logit"),
        offset = log(EFFORT_RATIO),
        weight = COMBINED_COUNT,
        data = x
      )
    
    bb3 <- 
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        dispformula = ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
        family = betabinomial(link = "logit"),
        offset = log(EFFORT_RATIO),
        weight = COMBINED_COUNT,
        data = x
      )
    
    binomial1 <-
      glmmTMB::glmmTMB(
        formula = PROP_15 ~ 1,
        family = binomial(link = "logit"),
        weight = COMBINED_COUNT,
        offset = log(EFFORT_RATIO),
        data = x
      )
    
    model_list <- 
      list(
        bb1 = bb1,
        bb2 = bb2,
        bb3 = bb3,
        binomial1 = binomial1
      )
    
    aic_table <- 
      make_aic_table(
        model_list 
      )
    
    loocv_table <-
      run_loocv(model_list = model_list, dat = dat, mapper = ccr_mapper)
    
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



fit_count_models <- 
  function(x, bootstrap_samples = NULL) {
    
    poisson1 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000),
        family = poisson(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        data = x
      )
    
    poisson2 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000) + I((CPUE_NO_KM2_15/1000)^2),
        family = poisson(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        data = x
      )
    
    nbin1 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000),
        family = glmmTMB::nbinom2(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        data = x
      )
    
    nbin2 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000),
        family = glmmTMB::nbinom2(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        disp = ~ LOG_CPUE_NO_KM2_15,
        data = x
      )
    
    nbin3 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000),
        family = glmmTMB::nbinom2(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        disp = ~ LOG_CPUE_NO_KM2_15 + I(LOG_CPUE_NO_KM2_15^2),
        data = x
      )
    
    nbin4 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000) + I((CPUE_NO_KM2_15/1000)^2),
        family = glmmTMB::nbinom2(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        data = x
      )
    
    nbin5 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000) + I((CPUE_NO_KM2_15/1000)^2),
        family = glmmTMB::nbinom2(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
        disp = ~ LOG_CPUE_NO_KM2_15,
        data = x
      )
    
    nbin6 <-
      glmmTMB::glmmTMB(
        formula = COUNT_30 ~ I(CPUE_NO_KM2_15/1000) + I((CPUE_NO_KM2_15/1000)^2),
        family = glmmTMB::nbinom2(link = "log"),
        offset = log(AREA_SWEPT_KM2_30),
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
      run_loocv(model_list = model_list, dat = dat, mapper = count_mapper)
    
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

# Mapping functions for converting model predictions to estimation scale

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

# Binomial and beta-binomial
ccr_mapper <- function(model, test_row, se_fit = FALSE) {
  
  if(se_fit) {
    # Two-fold CV
    inv_link  <- family(model)$linkinv
    
    p         <- predict(model, newdata = test_row, type = "link", se.fit = TRUE)
    p_lwr     <- p$fit - 2* p$se.fit
    p_upr     <- p$fit + 2* p$se.fit 
    
    ratio     <- (1 - p$fit) / p$fit
    ratio_lwr <- (1 - p_lwr) / p_lwr
    ratio_upr <- (1 - p_upr) / p_upr
    
    ratio     <- inv_link(ratio)
    ratio_lwr <- inv_link(ratio_lwr)
    ratio_upr <- inv_link(ratio_upr)
    
    output <- 
      data.frame(
        fit =     test_row$CPUE_NO_KM2_15 / ratio,
        fit_lwr = test_row$CPUE_NO_KM2_15 / ratio_lwr,
        fit_upr = test_row$CPUE_NO_KM2_15 / ratio_upr
      )
    
  } else {
    # LOOCV; Get predicted proportion (inv-logit scale)
    p <- predict(model, newdata = test_row, type = "response")
    ratio <- (1 - p) / p
    output <- test_row$CPUE_NO_KM2_15 / ratio
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
        fit_lwr = pred_count / test_row$AREA_SWEPT_KM2_30,
        fit_upr = pred_count / test_row$AREA_SWEPT_KM2_30
      )
    
  } else {
    # LOOCV; Get predicted counts
    pred_count <- predict(model, newdata = test_row, type = "response")
    output <- pred_count / test_row$AREA_SWEPT_KM2_30
  }

  
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



twofold_cv <- 
  function(model_list, validation_dat, mapper) {
    
    model_list = ols_results$models
    validation_dat = dat_oos
    mapper = ols_mapper
    
    obs_cpue    <- validation_dat$CPUE_NO_KM2_30
    obs_effort  <- validation_dat$AREA_SWEPT_KM2_30
    obs_count   <- validation_dat$COUNT_30
    
    results_list <- lapply(names(model_list), function(m_name) {
      mod <- model_list[[m_name]]
      
      preds <- mapper(
        model = mod,
        test_row = validation_dat,
        se_fit = TRUE
      )
      
      # Root mean square error
      output <- 
        data.frame(
          model_name = m_name,
          rmse       = sqrt(mean((preds$fit - obs_cpue)^2)),
          rmse_lci   = sqrt(mean((preds$fit_lwr - obs_cpue)^2)),
          rmse_uci   = sqrt(mean((preds$fit_upr - obs_cpue)^2)),
          tpe        = 100 * (sum(preds$fit * obs_effort) - sum(obs_count)) / sum(obs_count),
          tpe_lci    = 100 * (sum(preds$fit_lwr * obs_effort) - sum(obs_count)) / sum(obs_count),
          tpe_uci    = 100 * (sum(preds$fit_upr * obs_effort) - sum(obs_count)) / sum(obs_count)
        )
      
      output
    })
    
    results <- do.call(what = rbind, args = results_list)
    results$best <- results$rmse == min(results$rmse)
    
    return(results)
    
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

dharma_plots <- 
  function(model_list, common_name, subset = NULL, nsim = 1000, save_dir) {
    
    model_names    <- names(model_list)
    results        <- vector(mode = "list", length(model_list))
    names(results) <- model_names
    
    for(ii in seq_along(model_list)) {
      
      results[[ii]] <- DHARMa::simulateResiduals(model_list[[ii]], n = 1000)
      
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


# Fitting Somerton et al. (2002) 1998 data ---------------------------------------------------------

fpc_table_1998 <- data.frame()
loocv_table_1998 <- data.frame()
aic_table_1998 <- data.frame()
oos_table_1998 <- data.frame()

# Snow crab

# Set species
species_code <- 68580
common_name <- "snow crab"

# species_code <- 69322
# common_name <- "red king crab"

# species_code <- 68560
# common_name <- "Tanner crab"

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

# Fit models, check convergence
ols_results <- fit_ols_models(x = dat)
ccr_results <- fit_ccr_models(x = dat)
lognormal_results <- fit_lognormal(x = dat)
count_results <- fit_count_models(x = dat)

# Twofold CV
twofold_cv(model_list = ccr_results$models, validation_dat = dat_oos, mapper = ccr_mapper)
twofold_cv(model_list = lognormal_results$models, validation_dat = dat_oos, mapper = ratio_mapper)
twofold_cv(model_list = count_results$models, validation_dat = dat_oos, mapper = count_mapper)
twofold_cv(model_list = ols_results$models, validation_dat = dat_oos, mapper = ols_mapper)


# Check DHARMa residuals
dharma_plots(
  model_list = c(
    ols_results$models,
    ccr_results$models,
    lognormal_results$models,
    count_results$models
  ),
  common_name = common_name,
  nsim = 1000,
  subset = 1998,
  save_dir = dharma_dir
)

# Check OLS results
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

# Append cross validation results
loocv_table_1998 <- 
  dplyr::bind_rows(
    ols_results$loocv_table,
    ccr_results$loocv_table,
    lognormal_results$loocv_table,
    count_results$loocv_table
  ) |>
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(loocv_table_1998)


# Make AIC table that includes convergence checks
aic_table_1998 <-
  dplyr::bind_rows(
    ols_results$aic_table,
    ccr_results$aic_table,
    lognormal_results$aic_table,
    count_results$aic_table
  ) |>
  dplyr::bind_rows(aic_table_1998)

# Out-of-sample validation 
oos_table_1998 <- 
  indpendent_fpc_validation(
    fpc = fpc_1998,
    validation_dat = dat_oos
  ) |> 
  dplyr::arrange(rmse) |>
  dplyr::mutate(common_name = common_name) |>
  dplyr::bind_rows(oos_table_1998)
  


# Plot cross-year validation results

p_rmse_1998 <- 
  ggplot() +
  geom_point(
    data = oos_table_1998,
    mapping = aes(x = method, y = rmse),
    size = rel(2.2)) +
  geom_errorbar(
    data = oos_table_1998,
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
    data = oos_table_1998,
    mapping = aes(x = common_name, y = tpe, color = method), position = position_dodge(width = 0.5),
    size = rel(2.2)) +
  geom_errorbar(
    data = oos_table_1998,
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




# Fitting to the full data set ---------------------------------------------------------------------


