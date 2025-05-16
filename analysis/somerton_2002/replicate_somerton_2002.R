# Attempt to replicate Somerton et al. (2002)

# Somerton, D.A., Otto, R.S., Syrjala, S.E., 2002. Can changes in tow duration on bottom trawl surveys lead to changes in CPUE and mean size? Fish. Res. 55, 63–70. https://doi.org/10.1016/S0165-7836(01)00293-4

library(sratio)

somerton_crab <- readRDS(
  file = here::here("analysis", "somerton_2002", "data", "somerton_crab.rds")
  )

somerton_hauls_1998 <- readRDS(
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
  dplyr::mutate(CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
                COMBINED_COUNT = COUNT_30 + COUNT_15)


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

# Function to extract model intercept, variance, and bias-corrected ratio from Somerton's log-ratio models
bias_correction <- function(mod) {
  
  ratio_variance <- summary(mod)$sigma^2
  ratio_coef <- coef(mod)["(Intercept)"]
  ratio_estimate <- exp(ratio_coef + 0.5 * ratio_variance)
  
  output <- data.frame("ratio_intercept" = unname(ratio_coef),
    "ratio_variance" = unname(ratio_variance),
    "ratio_estimate" = unname(ratio_estimate))
  
  return(output)
  
}

# The best model for each species is an intercept-only model. Therefore, we use the intercept-only model for all three species.

# Estimated ratio without bias correction
raw_fit <- 
  data.frame(
    SPECIES_CODE = c(69322, 68580, 68560),
    ratio_estimate = exp(c(coef(mod_cpue_rkc_4)["(Intercept)"],
                           coef(mod_cpue_sc_4)["(Intercept)"],
                           coef(mod_cpue_tc_4)["(Intercept)"]
    )
    )
  )

# Estimated ratio with bias correction
somerton_ratios <- 
  cbind(
  data.frame(SPECIES_CODE = c(69322, 68580, 68560)),
  rbind(
    bias_correction(mod_cpue_rkc_4),
    bias_correction(mod_cpue_sc_4),
    bias_correction(mod_cpue_tc_4)
  )
)

# Make 1:1 lines
one_to_one <- 
  data.frame(SPECIES_CODE = c(69322, 68580, 68560),
             slope = 1,
             intercept = 0)

p1 <- 
  ggplot() +
  geom_abline(
    data = one_to_one,
    mapping = aes(color = "1:1 line",
                  linetype = "1:1 line",
                  slope = slope,
                  intercept = intercept)
  ) +
  geom_abline(
    data = raw_fit,
    mapping = 
      aes(
        slope = ratio_estimate,
        intercept = 0,
        color = "Ratio",
        linetype = "Ratio"
      )
  ) +
geom_abline(
  data = somerton_ratios, 
  mapping = 
    aes(
    slope = ratio_estimate, 
    intercept = 0,
    color = "Bias-corrected",
    linetype = "Bias-corrected"
  )
) +
  geom_point(
    data = cpue,
    mapping = aes(x = CPUE_NO_KM2_30, y = CPUE_NO_KM2_15),
    size = rel(0.6)
  ) +
  geom_rug(
    data = cpue,
    mapping = aes(x = CPUE_NO_KM2_30, y = CPUE_NO_KM2_15),
    sides = "b"
  ) +
  scale_color_manual(name = NULL, 
                     values = c("Bias-corrected" = "red",
                                             "Ratio" = "blue",
                                             "1:1 line" = "black")
                     ) +
  scale_linetype_manual(name = NULL, 
                        values = c("Bias-corrected" = 1,
                                                "Ratio" = 1,
                                                "1:1 line" = 2)
                        ) +
  scale_x_log10(name = expression(CPUE[30]~('#'%.%km^-2))) +
  scale_y_log10(name = expression(CPUE[15]~('#'%.%km^-2))) +
  facet_wrap(~sratio::species_code_label(x = SPECIES_CODE, type = "common_name"), scales = "free") +
  theme_bw() +
  theme(legend.title = element_blank())


png(filename = here::here("analysis", "somerton_2002", "plots", "compare_somerton.png"),
    width = 8, height = 3, units = "in", res = 300)
print(p1)
dev.off()
