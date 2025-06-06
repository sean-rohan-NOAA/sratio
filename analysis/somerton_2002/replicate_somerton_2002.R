library(sratio)

# Attempt to replicate Somerton et al. (2002)

# Somerton, D.A., Otto, R.S., Syrjala, S.E., 2002. Can changes in tow duration on bottom trawl surveys lead to changes in CPUE and mean size? Fish. Res. 55, 63–70. https://doi.org/10.1016/S0165-7836(01)00293-4

# Load data from Somerton's archived .txt data tables? If FALSE, uses filtered data that are available in CRABBASE. 
use_original_cpue_files = TRUE

somerton_crab <- 
  readRDS(
    file = here::here("analysis", "somerton_2002", "data", "somerton_crab.rds")
  )

somerton_hauls_1998 <- 
  readRDS(
    file = here::here("analysis", "somerton_2002", "data", "somerton_hauls_1998.rds")
  )

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

if(use_original_cpue_files) {
  
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
      CPUE_NO_KM2_30 = COUNT_30 / AREA_SWEPT_KM2_30
    ) |># Remove pairs with zero CPUE
    dplyr::mutate(
      LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
      LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
      CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
      CPUE_RATIO = CPUE_NO_KM2_15/CPUE_NO_KM2_30,
      COMBINED_COUNT = COUNT_30 + COUNT_15,
      common_name = sratio::species_code_label(SPECIES_CODE, type = "common_name")) |>
    dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0)
  
} else{
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
    dplyr::mutate(
      LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
      LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
      CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
      CPUE_RATIO = CPUE_NO_KM2_15/CPUE_NO_KM2_30,
      COMBINED_COUNT = COUNT_30 + COUNT_15,
      common_name = sratio::species_code_label(SPECIES_CODE, type = "common_name"))
}



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

par(mfrow = c(2,2))
plot(mod_cpue_tc_4)
plot(mod_cpue_rkc_4)
plot(mod_cpue_sc_4)

# Alternative: Zero-intercept linear regression between log-ratios
# Linear regression models might be suitable for log transformation
lm_rkc <- lm(LOG_CPUE_NO_KM2_15 ~ LOG_CPUE_NO_KM2_30 + 0, data = cpue[cpue$SPECIES_CODE == 69322, ])
lm_tc <- lm(LOG_CPUE_NO_KM2_15 ~ LOG_CPUE_NO_KM2_30 + 0, data = cpue[cpue$SPECIES_CODE == 68560, ])
lm_sc <- lm(LOG_CPUE_NO_KM2_15 ~ LOG_CPUE_NO_KM2_30 + 0, data = cpue[cpue$SPECIES_CODE == 68580, ])

plot(lm_rkc)
plot(lm_tc)
plot(lm_sc)

# Function to extract model intercept, variance, and bias-corrected ratio from log-ratio models
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

# Estimated ratio without bias correction using fixed effects ANOVA
ratios <- 
  cbind(
    data.frame(
      SPECIES_CODE = c(69322, 68580, 68560),
      common_name = sratio::species_code_label(x = c(69322, 68580, 68560), type = "common_name")),
    rbind(
      miller_bias_correct(mod_cpue_rkc_4), 
      miller_bias_correct(mod_cpue_sc_4),
      miller_bias_correct(mod_cpue_tc_4)
    )
  )

# Estimated ratio using a zero-intercept linear regression

fit_rkc <- miller_bias_correct(lm_rkc)
fit_rkc$fit$SPECIES_CODE <- 69322
fit_sc <- miller_bias_correct(lm_sc)
fit_sc$fit$SPECIES_CODE <- 68580
fit_tc <- miller_bias_correct(lm_tc)
fit_tc$fit$SPECIES_CODE <- 68560

zeroint_fit <-
  rbind(fit_rkc$fit,
        fit_sc$fit,
        fit_tc$fit)

zeroint_fit$common_name <- 
  sratio::species_code_label(
    x = zeroint_fit$SPECIES_CODE, 
    type = "common_name"
  )


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
