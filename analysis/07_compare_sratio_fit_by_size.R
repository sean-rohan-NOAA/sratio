library(sratio)

species_codes <- c(21740, 21720, 10110, 10112, 10115, 10130, 10210, 10261, 10285, 471, 68560, 68580, 69322)

ii <- 1

size_30 <- readRDS(here::here("output", "catch_at_length_1530.rds")) |>
  dplyr::filter(TREATMENT == 30, SPECIES_CODE == sp_codes[ii]) |>
  dplyr::ungroup()

size_15 <- readRDS(here::here("output", "catch_at_length_1530.rds")) |>
  dplyr::filter(TREATMENT == 15, SPECIES_CODE == sp_codes[ii]) |>
  dplyr::ungroup()

bootstrap_results_path <- list.files(here::here("output"), 
                                     recursive = TRUE, 
                                     pattern = "sratio_bootstrap_results_", 
                                     full.names = TRUE)

bootstrap_results <- readRDS(bootstrap_results_path[grep(pattern = sp_codes[ii], bootstrap_results_path)]) |>
  dplyr::group_by(SIZE_BIN) |>
  dplyr::summarise(s12 = mean(s12))

size_15 <- dplyr::select(size_30, MATCHUP, AREA_SWEPT_KM_30 = AREA_SWEPT_KM2) |> 
  unique() |>
  dplyr::inner_join(size_15) |>
  dplyr::inner_join(bootstrap_results)

# Predict catch-at-length, adjusting for area swept and selectivity
size_15$PREDICTED_FREQUENCY <- size_15$FREQ_EXPANDED * size_15$AREA_SWEPT_KM_30/size_15$AREA_SWEPT_KM2 * size_15$s12

# Predict catch-at-length adjusting for area swept but not selectivity
size_15$PREDICTED_FREQUENCY_NO_ADJ <- size_15$FREQ_EXPANDED * size_15$AREA_SWEPT_KM_30/size_15$AREA_SWEPT_KM2

compare_fit <- dplyr::select(size_15, MATCHUP, PREDICTED_FREQUENCY, PREDICTED_FREQUENCY_NO_ADJ, 
                         SIZE_BIN, SPECIES_CODE) |>
  dplyr::full_join(size_30)  |>
  dplyr::mutate(PREDICTED_FREQUENCY = if_else(is.na(PREDICTED_FREQUENCY), 0, PREDICTED_FREQUENCY),
                PREDICTED_FREQUENCY_NO_ADJ = if_else(is.na(PREDICTED_FREQUENCY_NO_ADJ), 0, PREDICTED_FREQUENCY_NO_ADJ),
                FREQ_EXPANDED = if_else(is.na(FREQ_EXPANDED), 0, FREQ_EXPANDED))

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


### Finish performance metrics by size bin


perf_bias
perf_mre
perf_mae
perf_r2

# Make plots
cowplot::plot_grid(p_fit_sratio_adjustment, p_area_only)




