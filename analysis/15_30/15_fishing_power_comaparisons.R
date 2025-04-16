# Catchability models for each taxa ----

library(sratio)
library(fishmethods)

seed <- 1729730909

# Load built-in data sets
catch_df <- sratio::data_1530$catch

haul_df <- sratio::data_1530$haul

sp_codes <- sort(unique(catch_df$SPECIES_CODE))


cpue_dat <- catch_df |>
  dplyr::inner_join(haul_df) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
                CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
  dplyr::mutate(LOG10_CPUE_30 = log10(CPUE_30+1),
                LOG10_CPUE_15 = log10(CPUE_15+1))

sp_codes <- sort(unique(catch_df$SPECIES_CODE))

dir.create(here::here("analysis", "15_30",  
                      "plots", "total_cpue_fit"),
           showWarnings = FALSE)


# Bayesian linear regression for log10(CPUE[30])~log10(CPUE[15]) ----

# May need to run species by species to avoid crashing
for(ii in 1:length(sp_codes)) {
  
  print(ii)
  
  sel_spp <- dplyr::filter(cpue_dat, SPECIES_CODE == sp_codes[ii])
  
  mod <- brms::brm(
    formula = LOG10_CPUE_30 ~ LOG10_CPUE_15 1, 
    data = sel_spp, 
    iter = 7000, 
    chains = 4,
    thin = 5,
    warmup = 2000
  )
  
  posterior_df <- brms::as_draws_df(mod)
  
  posterior_df$SPECIES_CODE <- sp_codes[ii]
  
  pp_check(mod, type = "dens_overlay", ndraws = 100)
  
  loo_check <- brms::loo(mod, moment_match = TRUE)
  
  loo::pareto_k_values(loo_check)
  
  saveRDS(object = posterior_df, 
          file = here::here("analysis", "15_30", "output", 
                            sp_codes[ii], 
                            paste0("cpue_brms_zeroint_posterior_", sp_codes[ii], ".rds")))
  
  pp_check(mod)
  
  rm(mod, posterior_df)
  
}

# test <-  lm(formula = I(CPUE_30/CPUE_15) ~ 1,
#    data = sel_spp)
# plot(test)

# hist((sel_spp$WEIGHT_30/sel_spp$AREA_SWEPT_KM2_30)/(sel_spp$WEIGHT_15/sel_spp$AREA_SWEPT_KM2_15))

zeroint_paths <- list.files(path = here::here("analysis", "15_30", "output"), 
                            pattern = "cpue_brms_zeroint_posterior", 
                            recursive = TRUE, 
                            full.names = TRUE)


cpue_fit <- data.frame()

for(jj in 1:length(zeroint_paths)) {
  
  cpue_fit <- dplyr::bind_rows(cpue_fit, readRDS(zeroint_paths[jj]))
  
}

slope_quantiles <- cpue_fit |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::summarise(qmin = min(b_LOG10_CPUE_15),
                   q025 = quantile(b_LOG10_CPUE_15, 0.025),
                   q250 = quantile(b_LOG10_CPUE_15, 0.25),
                   q500 = quantile(b_LOG10_CPUE_15, 0.5),
                   q750 = quantile(b_LOG10_CPUE_15, 0.75),
                   q975 = quantile(b_LOG10_CPUE_15, 0.975),
                   qmax = max(b_LOG10_CPUE_15),
                   type = "Slope")

saveRDS(object = slope_quantiles,
        file = here::here("analysis", "15_30", "output", "cpue_slope_posterior_quantiles.rds"))

# Plot results
png(here::here("analysis", "15_30", "plots", "total_cpue_fit", "cpue_model_density_plot.png"), 
    width = 169, height = 120, res = 300, units = "mm")
print(
ggplot() +
  geom_density(data = cpue_fit,
               mapping = aes(x = b_LOG10_CPUE_15)) +
  geom_vline(xintercept = 1, linetype = 2) +
  facet_wrap(~sratio:::species_code_label(SPECIES_CODE, type = "common_name", make_factor = TRUE), 
             scales = "free") +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = expression('Slope'~log[10](CPUE[30])*'~'*log[10](CPUE[15]))) +
  theme_bw() +
  theme(strip.background = element_blank())
)
dev.off()

png(here::here("analysis", "15_30", "plots", "total_cpue_fit", "cpue_model_violin_plot.png"), 
    width = 80, height = 100, res = 300, units = "mm")
print(
ggplot() +
  geom_violin(data = cpue_fit,
              mapping = aes(x = b_LOG10_CPUE_15, 
                            y = sratio:::species_code_label(SPECIES_CODE, 
                                                            type = "common_name", 
                                                            make_factor = TRUE)),
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(name = expression('Slope'~log[10](CPUE[30])*'~'*log[10](CPUE[15]))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9))
)
dev.off()

png(here::here("analysis", "15_30", "plots", "total_cpue_fit", "cpue_model_violin_plot.png"), 
    width = 80, height = 100, res = 300, units = "mm")
print(
  ggplot() +
    geom_violin(data = cpue_fit,
                mapping = aes(x = b_LOG10_CPUE_15, 
                              y = sratio:::species_code_label(SPECIES_CODE, 
                                                              type = "common_name", 
                                                              make_factor = TRUE)),
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(name = expression('Slope'~log[10](CPUE[30])*'~'*log[10](CPUE[15]))) +
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 9))
)
dev.off()


png(here::here("analysis", "15_30", "plots", "total_cpue_fit", "cpue_model_boxplot.png"), 
    width = 80, height = 60, res = 300, units = "mm")
print(
  ggplot() +
  geom_path(data = dplyr::group_by(cpue_fit, SPECIES_CODE) |>
              dplyr::summarise(q025 = quantile(b_LOG10_CPUE_15, 0.025),
                               q975 = quantile(b_LOG10_CPUE_15, 0.975)) |>
              tidyr::pivot_longer(cols = c(q025, q975)),
            mapping = aes(x = value, 
                          y = sratio:::species_code_label(SPECIES_CODE, type = "common_name", make_factor = TRUE))) +
  geom_path(data = dplyr::group_by(cpue_fit, SPECIES_CODE) |>
               dplyr::summarise(q250 = quantile(b_LOG10_CPUE_15, 0.25),
                                q750 = quantile(b_LOG10_CPUE_15, 0.75)) |>
              tidyr::pivot_longer(cols = c(q250, q750)),
             mapping = aes(x = value, 
                           y = sratio:::species_code_label(SPECIES_CODE, type = "common_name", make_factor = TRUE)),
            linewidth = 1.5) +
  geom_point(data = dplyr::group_by(cpue_fit, SPECIES_CODE) |>
               dplyr::summarise(median = median(b_LOG10_CPUE_15)),
             mapping = aes(x = median, 
                           y = sratio:::species_code_label(SPECIES_CODE, type = "common_name", make_factor = TRUE)), 
             shape = 21, fill = "white") +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_continuous(name = expression('Slope'~log[10](CPUE[30])*'~'*log[10](CPUE[15]))) +
  scale_y_discrete(limits = rev) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        axis.text = element_text(size = 8))
)
dev.off()


png(here::here("analysis", "15_30", "plots", "total_cpue_fit", "cpue_log_model_scatterplot.png"), 
    width = 169, height = 120, res = 300, units = "mm")
print(
ggplot(data = cpue_dat, 
       mapping = aes(x = log10(CPUE_15+1e-3), 
                     y = log10(CPUE_30+1e-3))) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ x + 0) +
  scale_x_continuous(name = expression(log[10](CPUE[15]))) +
  scale_y_continuous(name = expression(log[10](CPUE[30]))) +
  facet_wrap(~sratio:::species_code_label(SPECIES_CODE, type = "common_name", make_factor = TRUE), scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())
)
dev.off()

# CPUE ratio estimators from Wilderbuer (1998) ----
# Wilderbuer, T. K., R. F. Kappenman and D. R. Gunderson. 1998. Analysis of fishing power correction factor estimates from a trawl comparison experiment. North American Journal of Fisheries Management 18:11-18.

ratio_estimators <- 
  cpue_dat |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::summarise(
    ratio_rba = fishmethods::fpc(
      cpue1 = CPUE_30,
      cpue2 = CPUE_15,
      method = 2
    ),
    ratio_mm = fishmethods::fpc(
      cpue1 = CPUE_30,
      cpue2 = CPUE_15,
      method = 3
    ),
    ratio_kappenman = 
      fishmethods::fpc(
        cpue1 = CPUE_30,
        cpue2 = CPUE_15,
        method = 4
      )
  )


ratio_kappenman = 
  fishmethods::fpc(
    cpue1 = CPUE_30,
    cpue2 = CPUE_15,
    boot_type = "unpaired",
    kapp_zeros = "ind",
    method = 4
  )

rbind(
  cbind(ratio_estimators["SPECIES_CODE"], ratio_estimators[["ratio_rba"]]),
  cbind(ratio_estimators["SPECIES_CODE"], ratio_estimators[["ratio_mm"]]),
  cbind(ratio_estimators["SPECIES_CODE"], ratio_estimators[["ratio_kappenman"]])
)


ratio_randomized_block_anova = fishmethods::fpc(
  cpue1 = cpue_dat$CPUE_30,
  cpue2 = cpue_dat$CPUE_15,
  method = 2
)


# Calculate mean bias ----
# CIs estimated using haul-level bootstrap

# Load built-in data sets
catch_df <- sratio::data_1530$catch #|>
# dplyr::filter(CRUISE %in% use_cruises)

haul_df <- sratio::data_1530$haul #|>
# dplyr::filter(CRUISE %in% use_cruises)

sp_codes <- sort(unique(catch_df$SPECIES_CODE))

bias_table <- catch_df |>
  dplyr::inner_join(haul_df) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
                CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
  dplyr::mutate(ratio = CPUE_15/CPUE_30,
                log_error = log10(CPUE_30+1)-log10(CPUE_15+1),
                abs_error = abs(log10(CPUE_30+1)-log10(CPUE_15+1)),
                sq_error = (CPUE_30-CPUE_15)^2) |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::reframe(MEAN_RATIO = mean(ratio),
                   BIAS = 10^(mean(log_error)),
                   MAE = 10^mean(abs_error),
                   RMSE = sqrt(mean(sq_error))) |>
  dplyr::mutate(COMMON_NAME = sratio:::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::select(COMMON_NAME, BIAS, MEAN_RATIO, MAE, RMSE)

bias_samples <- catch_df |>
  dplyr::inner_join(haul_df) |>
  dplyr::select(MATCHUP, SPECIES_CODE, WEIGHT, MATCHUP, TREATMENT) |>
  dplyr::mutate(TREATMENT = paste0("WEIGHT_", TREATMENT)) |>
  tidyr::pivot_wider(names_from = TREATMENT, values_from = WEIGHT, values_fill = 0) |>
  dplyr::inner_join(readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
                      dplyr::select(SPECIES_CODE, MATCHUP, AREA_SWEPT_KM2_15, AREA_SWEPT_KM2_30) |>
                      unique() |>
                      dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)))) |>
  dplyr::mutate(CPUE_30 = WEIGHT_30/AREA_SWEPT_KM2_30,
                CPUE_15 = WEIGHT_15/AREA_SWEPT_KM2_15) |>
  dplyr::group_by(SPECIES_CODE) |>
  # Bootstrap
  dplyr::reframe(BIAS = bootstrap_mean_bias(CPUE_30, CPUE_15, n_samples = 10000, add_constant = 1, scale = "log10", seed = seed)) |>
  dplyr::mutate(COMMON_NAME = sratio:::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::select(SPECIES_CODE, BIAS)

bias_quantiles <- bias_samples |>
  dplyr::group_by(SPECIES_CODE) |>
  dplyr::summarise(q025 = quantile(BIAS, 0.025),
                   q250 = quantile(BIAS, 0.25),
                   q500 = quantile(BIAS, 0.5),
                   q750 = quantile(BIAS, 0.75),
                   q975 = quantile(BIAS, 0.975),
                   type = "Bias")

saveRDS(bias_quantiles, 
        file = here::here("analysis", "15_30", "output", "cpue_bias_bootstrap_quantiles.rds"))


lines <- c("CPUE comparison between 15 and 30 minute tows\n", "BIAS > 1 = 30 minutes higher\n\n\n")
cat(lines, file = here::here("analysis", "15_30", "plots", "total_cpue_fit", "bias_table.csv"))

write.table(bias_table, 
            file = here::here("analysis", "15_30", "plots", "total_cpue_fit", "bias_table.csv"), 
            append = TRUE, 
            row.names = FALSE, 
            sep = ",")
