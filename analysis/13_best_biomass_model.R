library(sratio)

sp_codes <- sort(unique(sratio::data_1530$catch$SPECIES_CODE))

treatments <- factor(c(15, 30))

n_cores <- 4

aic_df <- data.frame()

for(ii in 1:length(sp_codes)) {
  
  sel_dat <- merge(
    sratio::data_1530$catch[sratio::data_1530$catch$SPECIES_CODE == sp_codes[ii], ],
    sratio::data_1530$haul[c("VESSEL", "CRUISE", "HAUL", "TREATMENT", "MATCHUP", "AREA_SWEPT_KM2", "DURATION", "HAULJOIN")],
    all.y = TRUE
  )
  
  # Only use data from hauls where at least one matchup had a positive catch rate
  na_count <- table(is.na(sel_dat$NUMBER_FISH), 
                    sel_dat$MATCHUP)
  
  na_matchup <- as.numeric(names(which(na_count[2, ] > 1)))
  
  sel_dat <- sel_dat[!(sel_dat$MATCHUP %in% na_matchup), ]
  
  sel_dat$WEIGHT[is.na(sel_dat$WEIGHT)] <- 0
  sel_dat$SPECIES_CODE[is.na(sel_dat$SPECIES_CODE)] <- sp_codes[ii]
  sel_dat$NUMBER_FISH[is.na(sel_dat$NUMBER_FISH)] <- 0
  sel_dat$MATCHUP <- factor(sel_dat$MATCHUP)
    
  m0 <- mgcv::gam(formula = WEIGHT ~ s(MATCHUP, bs = "re") + offset(I(log(AREA_SWEPT_KM2))), 
                  data = sel_dat,
                  family = gaussian(link = "log")) 
  
  m0_gamma <- mgcv::gam(formula = I(WEIGHT+1e-5) ~ s(MATCHUP, bs = "re") + offset(I(log(AREA_SWEPT_KM2))), 
                        data = sel_dat,
                        family = Gamma(link = "log")) 
  
  m1 <- mgcv::gam(formula = WEIGHT ~ 0 + TREATMENT + s(MATCHUP, bs = "re") + offset(I(log(AREA_SWEPT_KM2))), 
                  data = sel_dat,
                  family = gaussian(link = "log"))
  
  m1_gamma <- mgcv::gam(formula = I(WEIGHT+1e-5) ~ 0 + TREATMENT + s(MATCHUP, bs = "re") + offset(I(log(AREA_SWEPT_KM2))), 
                        data = sel_dat,
                        family = Gamma(link = "log"))
  
  aic_table <- AIC(m0_gamma, m1_gamma)
  
  aic_table$FORMULA <- as.character(
    c(
      m0_gamma$call$formula, 
      m1_gamma$call$formula
    )
  )
  
  aic_table$SPECIES_CODE <- sp_codes[ii]
  aic_table$MODEL <- rownames(aic_table)
  
  aic_table <- aic_table[c("SPECIES_CODE", "MODEL", "AIC", "df", "FORMULA")]
  
  aic_table <- aic_table[order(aic_table$AIC), ]
  
  aic_lowest <- which.min(aic_table$AIC)
  
  aic_within_2 <- which(aic_table$AIC < min(aic_table$AIC) + 2)
  
  aic_better_k <- which(aic_table$df[aic_within_2] < aic_table$df[aic_lowest])
  
  if(length(aic_better_k) == 0) { 
    aic_best <- aic_lowest 
    } else {
    aic_best <- aic_better_k[aic_table$df[aic_better_k] == min(aic_table$df[aic_better_k])]
  }
         
  if(length(aic_best) > 1) {
    aic_best <- aic_best[which.min(aic_table$AIC[aic_best])][1]
  }            
  
  aic_table$BEST_MODEL <- 1:nrow(aic_table) == aic_best

  
  aic_df <- rbind(aic_table, aic_df)
}


aic_df[aic_df$BEST_MODEL, ]





cpue_dat <- sratio::data_1530$catch |>
  dplyr::inner_join(sratio::data_1530$haul) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("output", "n_by_treatment_1530.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
                CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
  dplyr::mutate(LOG10_CPUE_30 = log10(CPUE_30+0.001),
                LOG10_CPUE_15 = log10(CPUE_15+0.001))


sp_codes <- sort(unique(sratio::data_1530$catch$SPECIES_CODE))

for(ii in 1:length(sp_codes)) {
  
  sel_spp <- dplyr::filter(cpue_dat, SPECIES_CODE == sp_codes[ii])
  
  mod <- brms::brm(formula = LOG10_CPUE_15 ~ LOG10_CPUE_30 + 0, 
                   data = sel_spp, 
                   iter = 3000, 
                   chains = 4, 
                   warmup = 1000)
  
  posterior_df <- brms::as_draws_df(mod)
  
  posterior_df$SPECIES_CODE <- sp_codes[ii]
  
}



ggplot() +
  geom_density(data = posterior_df,
               mapping = aes(x = b_LOG10_CPUE_30))



ggplot(data = test, mapping = aes(x = CPUE_30, y = CPUE_15)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_x_continuous(name = expression(CPUE[30]~(kg/km^2))) +
  scale_y_continuous(name = expression(CPUE[15]~(kg/km^2))) +
  facet_wrap(~sratio:::species_code_label(SPECIES_CODE, type = "common_name"), scales = "free") +
  theme_bw()


ggplot(data = test, mapping = aes(x = log10(CPUE_30+1e-3), y = log10(CPUE_15+1e-3))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_x_continuous(name = expression(log[10](CPUE[30])~(kg/km^2))) +
  scale_y_continuous(name = expression(log[10](CPUE[15])~(kg/km^2))) +
  facet_wrap(~sratio:::species_code_label(SPECIES_CODE, type = "common_name"), scales = "free") +
  theme_bw()




