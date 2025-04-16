# Attempt to replicate Somerton et al. (2002)

# Somerton, D.A., Otto, R.S., Syrjala, S.E., 2002. Can changes in tow duration on bottom trawl surveys lead to changes in CPUE and mean size? Fish. Res. 55, 63–70. https://doi.org/10.1016/S0165-7836(01)00293-4

library(sratio)

crab_species_codes <- c(68560, 68580, 69322)

channel <- sratio::get_connected()

somerton_hauls_1998 <- 
  RODBC::sqlQuery(
    channel = channel,
    query = 
      paste0(
        "SELECT hauljoin, vessel, cruise, haul, duration, distance_fished, net_width 
      FROM racebase.haul
      WHERE
        cruise = 199801
        AND hauljoin IN (", paste(unique(sratio::crab_size_1995_1998$HAULJOIN), collapse = ", "),
        ")"
      ) 
  ) |>
  dplyr::mutate(
    AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
    TREATMENT = 
      factor(
        plyr::round_any(DURATION, 0.25), 
        levels = c(0.25, 0.5), 
        labels = c(15, 30)
      )
  ) |>
  dplyr::inner_join(
    dplyr::select(
      sratio::otto_key_1998, # Use tow pairs from Somerton et al. (2002)
      VESSEL, CRUISE, HAUL, TOW_PAIR, TOW_BLOCK
    ),
    by = c("VESSEL", "CRUISE", "HAUL")
  ) |>
  dplyr::arrange(TOW_BLOCK)

somerton_crab <- 
  somerton_hauls_1998 |>
  dplyr::select(HAULJOIN, VESSEL, CRUISE, HAUL, AREA_SWEPT_KM2, DURATION, TREATMENT, TOW_PAIR, TOW_BLOCK) |>
  dplyr::inner_join(
    dplyr::filter(sratio::crab_size_1995_1998),
    by = c("HAULJOIN", "VESSEL", "CRUISE", "HAUL")
    ) |>
      dplyr::filter(SPECIES_CODE %in% crab_species_codes,  # RKC, snow, Tanner from 1998
                    CRUISE == 199801) |>
      dplyr::mutate(SPECIES_CODE_SEX  = as.numeric(paste0(SPECIES_CODE, SEX)),  # Combine sex and species code
                    VESSEL = factor(VESSEL),
                    SEX = factor(SEX)) # Make vessel a factor for models

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


  test <- cpue[cpue$SPECIES_CODE == 69322, ]

  mod_cpue_rkc <- 
  lm(
    formula = I(log(CPUE_NO_KM2_15/CPUE_NO_KM2_30)) ~ SEX + VESSEL,
    data = test
  )
  
  cov(test$CPUE_NO_KM2_15, test$CPUE_NO_KM2_30)

summary(mod_cpue_rkc)
anova(mod_cpue_rkc)

mod_cpue_sc <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX + VESSEL,
    data = cpue[cpue$SPECIES_CODE == 68580, ]
  )

summary(mod_cpue_sc)
anova(mod_cpue_sc)

mod_cpue_tc <- 
  lm(
    formula = CPUE_LOG_RATIO ~ SEX + VESSEL,
    data = cpue[cpue$SPECIES_CODE == 68560, ]
  )

summary(mod_cpue_tc)
anova(mod_cpue_tc)


test <- glm(
  formula = cbind(COUNT_15, COUNT_30) ~ 1 + offset(log(AREA_SWEPT_KM2_15/AREA_SWEPT_KM2_30)),
  family = binomial(),
  data = dplyr::filter(cpue, SPECIES_CODE == 68560, SEX == 2)
)

test$coefficients
sratio::inv_logit(test$coefficients)


predict(test, newdata = data.frame(AREA_SWEPT_KM2_15 = 0.024, AREA_SWEPT_KM2_30 = c(0.048)), type = "response")

# # What if sampling factors weren't used -----
# # THIS IS NOT THE CORRECT APPROACH
# 
# cpue_wrong <- 
#   somerton_crab |>
#   dplyr::group_by(VESSEL, CRUISE, HAUL, TOW_PAIR, TREATMENT, AREA_SWEPT_KM2, SPECIES_CODE, SEX) |>
#   dplyr::summarize(COUNT = n()) |>
#   dplyr::ungroup() |>
#   dplyr::mutate(CPUE_NO_KM2 = COUNT/AREA_SWEPT_KM2) |>
#   dplyr::select(-HAUL, -CRUISE) |>
#   tidyr::pivot_wider(names_from = "TREATMENT", 
#                      values_from = c("CPUE_NO_KM2", "COUNT", "AREA_SWEPT_KM2"),
#                      values_fill = 0) |>
#   dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0) |> # Remove pairs with zero CPUE
#   dplyr::mutate(CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30))
# 
# 
# mod_cpue_wrong_rkc <- 
#   lm(
#     formula = CPUE_LOG_RATIO ~ 1,
#     data = cpue_wrong[cpue_wrong$SPECIES_CODE == 69322, ]
#   )
# 
# summary(mod_cpue_wrong_rkc)
# anova(mod_cpue_wrong_rkc)
# 
# mod_cpue_wrong_sc <- 
#   lm(
#     formula = CPUE_LOG_RATIO ~ 1,
#     data = cpue_wrong[cpue_wrong$SPECIES_CODE == 68580, ]
#   )
# 
# summary(mod_cpue_wrong_sc)
# anova(mod_cpue_wrong_sc)
# 
# mod_cpue_wrong_tc <- 
#   lm(
#     formula = CPUE_LOG_RATIO ~ 1,
#     data = cpue[cpue$SPECIES_CODE == 68560, ]
#   )

# Function to extract model intercept, variance, and bias-corrected ratio from Somerton's log-ratio models
extract_bias_corrected_ratio <- function(mod) {
  
  ratio_variance <- summary(mod)$sigma^2
  ratio_coef <- coef(mod)["(Intercept)"]
  ratio_estimate <- exp(ratio_coef + 0.5 * ratio_variance)
  
  return(
    c("ratio_intercept" = unname(ratio_coef),
      "ratio_variance" = unname(ratio_variance),
      "ratio_estimate" = unname(ratio_estimate))
  )
  
}

extract_bias_corrected_ratio(mod_cpue_wrong_tc)
extract_bias_corrected_ratio(mod_cpue_tc)
extract_bias_corrected_ratio(mod_cpue_sc)
extract_bias_corrected_ratio(mod_cpue_rkc)



fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68580 & cpue$SEX == 1],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580 & cpue$SEX == 1],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68580 & cpue$SEX == 2],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68580 & cpue$SEX == 2],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 69322 & cpue$SEX == 1],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322 & cpue$SEX == 1],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 69322 & cpue$SEX == 2],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 69322 & cpue$SEX == 2],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68560 & cpue$SEX == 1],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560 & cpue$SEX == 1],
                 method = 4)

fishmethods::fpc(cpue1 = cpue$CPUE_NO_KM2_15[cpue$SPECIES_CODE == 68560 & cpue$SEX == 2],
                 cpue2 = cpue$CPUE_NO_KM2_30[cpue$SPECIES_CODE == 68560 & cpue$SEX == 2],
                 method = 4)





ggplot() +
  geom_point(data = cpue_wrong,
             mapping = aes(x = CPUE_NO_KM2_30, y = CPUE_NO_KM2_15)) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~SPECIES_CODE)


cowplot::plot_grid(
  ggplot() +
    geom_histogram(data = cpue,
                   mapping = aes(x = CPUE_NO_KM2_15/CPUE_NO_KM2_30)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~SPECIES_CODE),
  ggplot() +
    geom_histogram(data = cpue_wrong,
                   mapping = aes(x = CPUE_NO_KM2_15/CPUE_NO_KM2_30)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~SPECIES_CODE),
  nrow = 2
)


ggplot() +
  stat_ecdf(data = cpue,
                 mapping = aes(x = CPUE_NO_KM2_15/CPUE_NO_KM2_30)) +
  # scale_x_log10() +
  # scale_y_log10() +
  facet_wrap(~SPECIES_CODE)
