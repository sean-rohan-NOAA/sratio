library(sratio)

# Load built-in data sets
catch_df <- sratio::data_1530$catch #|>
  # dplyr::filter(CRUISE %in% use_cruises)

haul_df <- sratio::data_1530$haul #|>
  # dplyr::filter(CRUISE %in% use_cruises)

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

# May need to run species by species to avoid crashing --- 
for(ii in 1:length(sp_codes)) {
  
  print(ii)
  
  sel_spp <- dplyr::filter(cpue_dat, SPECIES_CODE == sp_codes[ii])
  
  mod <- brms::brm(formula = LOG10_CPUE_30 ~ LOG10_CPUE_15 + 0, 
                   data = sel_spp, 
                   iter = 5000, 
                   chains = 4,
                   warmup = 2000)
  
  posterior_df <- brms::as_draws_df(mod)
  
  posterior_df$SPECIES_CODE <- sp_codes[ii]
  
  saveRDS(object = posterior_df, 
          file = here::here("analysis", "15_30", "output", 
                            sp_codes[ii], 
                            paste0("cpue_brms_zeroint_posterior_", sp_codes[ii], ".rds")))
  
  rm(mod, posterior_df)
  
}

zeroint_paths <- list.files(path = here::here("analysis", "15_30", "output"), 
                            pattern = "cpue_brms_zeroint_posterior", 
                            recursive = TRUE, 
                            full.names = TRUE)


cpue_fit <- data.frame()

for(jj in 1:length(zeroint_paths)) {
  
  cpue_fit <- dplyr::bind_rows(cpue_fit, readRDS(zeroint_paths[jj]))
  
}

# Plot CPUE model results -----
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

# cpue_plot_dat <- cpue_fit |>
#   dplyr::group_by(SPECIES_CODE) |>
#   dplyr::summarise(qmin = min(b_LOG10_CPUE_15),
#                    q025 = quantile(b_LOG10_CPUE_15, 0.025),
#                    q250 = quantile(b_LOG10_CPUE_15, 0.25),
#                    q500 = quantile(b_LOG10_CPUE_15, 0.5),
#                    q750 = quantile(b_LOG10_CPUE_15, 0.75),
#                    q975 = quantile(b_LOG10_CPUE_15, 0.975),
#                    qmax = max(b_LOG10_CPUE_15))


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

