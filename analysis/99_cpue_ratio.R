# Model 


library(lme4)
library(ggplot2)
library(tidyr)

dat <- dplyr::inner_join(readRDS(file = here::here("data", "catch_1530.rds")),
                         readRDS(file = here::here("data", "hauls_1530.rds")), 
                         by = c("VESSEL", "CRUISE", "HAUL", "HAULJOIN", "MATCHUP")) |>
  dplyr::mutate(TREATMENT = factor(TREATMENT)) |>
  dplyr::mutate(CPUE_KGKM2 = WEIGHT/AREA_SWEPT_KM2,
                CPUE_NKM2  = NUMBER_FISH/AREA_SWEPT_KM2,
                YEAR = lubridate::year(START_TIME))

ratio_3015 <- dat |>
  dplyr::select(MATCHUP, CPUE_KGKM2, TREATMENT, SPECIES_CODE, YEAR) |>
  tidyr::pivot_wider(names_from = "TREATMENT", values_from = "CPUE_KGKM2") |>
  dplyr::mutate(RATIO_1530 = `15`/`30`)

ggplot() +
  geom_point(data = ratio_3015,
             mapping = aes(x = MATCHUP, 
                           y = RATIO_1530)) +
  geom_hline(yintercept = 1) +
  facet_wrap(~SPECIES_CODE, ncol = 2, scales = "free_y")

ggplot() +
  geom_density(data = ratio_3015,
               mapping = aes(x = log10(RATIO_1530))) +
  geom_vline(xintercept = 0) +
  facet_wrap(~SPECIES_CODE, ncol = 2, scales = "free_x")

ggplot() +
  geom_histogram(data = ratio_3015,
                 mapping = aes(x = RATIO_1530)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~SPECIES_CODE, ncol = 2, scales = "free_x") +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10))


dat |>
  dplyr::select(MATCHUP, CPUE_KGKM2, TREATMENT, SPECIES_CODE, YEAR) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = CPUE_KGKM2) |>
  ggplot() +
  geom_point(mapping = aes(x = `30`, y = `15`)) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~SPECIES_CODE, scales = "free")


fit_cpue_weight_models <- function(catch_weight, 
                                   log_offset = 0.001, 
                                   area_swept = NULL, 
                                   haul_duration = NULL, 
                                   treatment, 
                                   matchup, 
                                   family = gaussian()) {
  
  m1 <- m2 <- m3 <- m4 <- m5 <- m6 <- NULL

  input_data <- data.frame(catch_weight = catch_weight,
                           area_swept = area_swept,
                           haul_duration = haul_duration,
                           treatment = treatment,
                           matchup = matchup,
                           log_offset = log_offset,
                           dummy_var = 1)

  if(!is.null(haul_duration)) {
    m1 <- mgcv::gam(formula = log(catch_weight + log_offset) ~ treatment + 
                      s(matchup, bs = "re", by = dummy_var) + 0, 
                     offset = haul_duration,
                     family = family,
                     data = input_data)
    
    m3 <- mgcv::gam(formula = log(catch_weight + log_offset) ~ treatment + 
                      s(log(haul_duration), bs = "tp") + 
                      s(matchup, bs = "re", by = dummy_var) + 0, 
                     family = family,
                     data = input_data)
    
    m5 <- mgcv::gam(formula = log(catch_weight + log_offset) ~ treatment + 
                       s(log(haul_duration), bs = "tp") + 
                       s(matchup, bs = "re", by = dummy_var) + 0, 
                     offset = haul_duration,
                     family = family,
                     data = input_data)
  }

  if(!is.null(area_swept)) {
    m2 <- mgcv::gam(formula = log(catch_weight + log_offset) ~ treatment + 
                      s(matchup, bs = "re", by = dummy_var) + 0, 
                     offset = area_swept,
                     family = family,
                     data = input_data)
    
    m4 <- mgcv::gam(formula = log(catch_weight + log_offset) ~ treatment + 
                      s(log(area_swept), bs = "tp") + 
                      s(matchup, bs = "re", by = dummy_var) + 0, 
                     family = family,
                     data = input_data)
    
    m6 <- mgcv::gam(formula = log(catch_weight + log_offset) ~ treatment + 
                      s(log(area_swept), bs = "tp") + 
                      s(matchup, bs = "re", by = dummy_var) + 0, 
                     offset = area_swept,
                     family = family,
                     data = input_data)
  }
  
  out <- list(models = list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5, m6 = m6), 
              input_data = input_data)
  
  return(out)

  
}


input_species <- unique(dat$SPECIES_CODE)

for(ii in 1:length(input_species)) {
  
  sel_dat <- dplyr::filter(dat, SPECIES_CODE == input_species[ii])
  
  spp_models <- fit_cpue_weight_models(catch_weight = sel_dat$WEIGHT, 
                                       log_offset = 0.001, 
                                       area_swept = sel_dat$AREA_SWEPT_KM2, 
                                       haul_duration = sel_dat$DURATION,
                                       matchup = sel_dat$MATCHUP,
                                       treatment = sel_dat$TREATMENT)
  
  assign(x = paste0("cpue_models_", input_species[ii]),
         value = spp_models)
  
}


exp(coef(cpue_models_10110$models$m1$gam))
exp(coef(cpue_models_10130$models$m1$gam))
exp(coef(cpue_models_10210$models$m1$gam))
exp(coef(cpue_models_10261$models$m1$gam))
exp(coef(cpue_models_21720$models$m1$gam))
exp(coef(cpue_models_21740$models$m1$gam))
exp(coef(cpue_models_471$models$m1$gam))
exp(coef(cpue_models_68560$models$m1$gam))
exp(coef(cpue_models_68580$models$m1$gam))
exp(coef(cpue_models_69322$models$m1$gam))

exp(coef(cpue_models_10110$models$m2$gam))
exp(coef(cpue_models_10130$models$m2$gam))
exp(coef(cpue_models_10210$models$m2$gam))
exp(coef(cpue_models_10261$models$m2$gam))
exp(coef(cpue_models_21720$models$m2$gam))
exp(coef(cpue_models_21740$models$m2$gam))
exp(coef(cpue_models_471$models$m2$gam))
exp(coef(cpue_models_68560$models$m2$gam))
exp(coef(cpue_models_68580$models$m2$gam))
exp(coef(cpue_models_69322$models$m2$gam))

ggplot() +
  geom_point(data = dat,
             mapping = aes(x = MATCHUP, 
                           y = log10(CPUE_KGKM2), 
                           color = factor(TREATMENT))) +
  facet_wrap(~SPECIES_CODE, ncol = 2, scales = "free_y")


