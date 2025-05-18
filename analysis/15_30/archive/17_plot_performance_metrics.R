library(sratio)

performance_metrics <- dplyr::bind_rows(
  readRDS(here::here("analysis", "15_30", "output", "cpue_slope_posterior_quantiles.rds")),
  readRDS(here::here("analysis", "15_30", "output", "cpue_bias_bootstrap_quantiles.rds"))
)


performance_metrics$COMMON_NAME <- sratio:::species_code_label(performance_metrics$SPECIES_CODE, 
                                                               type = "common_name", 
                                                               make_factor = TRUE)

axis_ticks_labels <- data.frame(COMMON_NAME = unique(performance_metrics$COMMON_NAME)) |>
  dplyr::mutate(Y = rank(as.numeric(COMMON_NAME)))

intermediate_ticks <- data.frame(Y = seq(min(axis_ticks_labels$Y) - 0.5, max(axis_ticks_labels$Y) + 0.5, 1))

performance_metrics <- performance_metrics |>
  dplyr::inner_join(axis_ticks_labels) |>
  dplyr::mutate(Y = if_else(type == "Slope", Y + 0.25, Y - 0.25))

png(here::here("analysis", "15_30", "plots", "total_cpue_fit", "cpue_metric_boxplots.png"), 
    width = 80, height = 100, res = 300, units = "mm")
print(
  ggplot() +
    geom_hline(data = intermediate_ticks,
               mapping = aes(yintercept = Y),
               color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = 2, color = "red") +
    geom_segment(data = performance_metrics,
                 mapping = aes(x = q025,
                               xend = q975,
                               y = Y,
                               color = type)) +
    geom_segment(data = performance_metrics,
                 mapping = aes(x = q250,
                               xend = q750,
                               y = Y,
                               color = type),
                 linewidth = 1.5) +
    geom_point(data = performance_metrics,
               mapping = aes(x = q500, 
                             y = Y,
                             color = type), 
               shape = 21, fill = "white") +
    annotate(label = "15 higher",
             x = -Inf,
             y = Inf,
             geom = "text",
             hjust = -0.1,
             vjust = -0.5,
             color = "red") +
    annotate(label = "30 higher",
             x = Inf,
             y = Inf,
             geom = "text",
             hjust = 1.1,
             vjust = -0.5,
             color = "red") +
    scale_y_reverse(breaks = rev(axis_ticks_labels$Y), labels = rev(axis_ticks_labels$COMMON_NAME)) +
    scale_x_continuous(name = "Value") +
    scale_color_colorblind(name = "Metric") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 8),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom")
)
dev.off()

perf <- readRDS(file = here::here("analysis", "15_30", "output", 
                                  paste0("model_performance_", model_method, "_", sp_codes[ii], ".rds")))


compare_fit_long <- tidyr::pivot_longer(compare_fit, 
                                        cols = c("PREDICTED_FREQUENCY", "PREDICTED_FREQUENCY_NO_ADJ"))

p_fit_sratio_adjustment <- ggplot(data = compare_fit,
                                  mapping = aes(x = SIZE_BIN,
                                                y = PREDICTED_FREQUENCY-FREQ_EXPANDED)) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(name = expression("Selectivity ratio "~hat(N[30])-N[30])) +
  scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
  theme_bw()

p_area_only <- ggplot(data = compare_fit,
                      mapping = aes(x = SIZE_BIN,
                                    y = PREDICTED_FREQUENCY_NO_ADJ-FREQ_EXPANDED)) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(name = expression("No SR adjustment "~hat(N[30])-N[30])) +
  scale_x_continuous(name = sratio::species_code_label(x = sp_codes[ii])) +
  theme_bw()

p_fit_vs_obs <- ggplot() +
  geom_point(data = compare_fit,
             mapping = aes(x = PREDICTED_FREQUENCY,
                           y = FREQ_EXPANDED)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_x_continuous(name = expression("Selectivity ratio fit "~hat(N[30]))) + 
  scale_y_continuous(name = expression("Observed "~N[30])) +
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