
perf <- readRDS(file = here::here("analysis", "shelf_slope", "output", 
                                  paste0("model_performance_", model_method, "_", sp_codes[ii], ".rds")))


compare_fit_long <- tidyr::pivot_longer(compare_fit, 
                                        cols = c("PREDICTED_FREQUENCY", "PREDICTED_FREQUENCY_NO_ADJ"))

p_fit_sratio_adjustment <- ggplot(data = compare_fit,
                                  mapping = aes(x = SIZE_BIN,
                                                y = PREDICTED_FREQUENCY-FREQ_EXPANDED)) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(name = expression("Selectivity ratio "~hat(N[PNE])-N['83-112'])) +
  scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
  theme_bw()

p_area_only <- ggplot(data = compare_fit,
                      mapping = aes(x = SIZE_BIN,
                                    y = PREDICTED_FREQUENCY_NO_ADJ-FREQ_EXPANDED)) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(name = expression("No SR adjustment "~hat(N[PNE])-N['83-112'])) +
  scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
  theme_bw()

p_fit_vs_obs <- ggplot() +
  geom_point(data = compare_fit,
             mapping = aes(x = PREDICTED_FREQUENCY,
                           y = FREQ_EXPANDED)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_continuous(name = expression("Selectivity ratio fit "~hat(N[PNE]))) + 
  scale_y_continuous(name = expression("Observed "~N[PNE])) +
  facet_wrap(~SIZE_BIN, scale = "free") +
  theme_bw()

p_bias <- ggplot(data = compare_fit,
                 mapping = aes(x = SIZE_BIN,
                               y = 10^(log10(PREDICTED_FREQUENCY)-log10(FREQ_EXPANDED)))) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(name = "Bias") +
  scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
  theme_bw()

p_mre <- ggplot(data = compare_fit,
                mapping = aes(x = SIZE_BIN,
                              y = (PREDICTED_FREQUENCY-FREQ_EXPANDED)/FREQ_EXPANDED*100)) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(name = "MRE (%)") +
  scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
  theme_bw()

# Make plots
cowplot::plot_grid(p_fit_sratio_adjustment, p_area_only)
