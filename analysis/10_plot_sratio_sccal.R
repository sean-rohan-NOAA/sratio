# Plot catch-at-length model results
library(sratio)

sccal_bootstrap_path <- list.files(here::here("output"), 
                                     recursive = TRUE, 
                                     pattern = "sccal_model_bootstrap_results_", 
                                     full.names = TRUE)

sratio_bootstrap_path <- list.files(here::here("output"), 
                                   recursive = TRUE, 
                                   pattern = "sratio_bootstrap_results_", 
                                   full.names = TRUE)
sp_codes <- unique(c(as.numeric(gsub("[^0-9]", "", basename(sratio_bootstrap_path))),
                     as.numeric(gsub("[^0-9]", "", basename(sccal_bootstrap_path))))
                   )
as.numeric(gsub("[^0-9]", "", basename(sccal_bootstrap_path)))

obs_ratio <- readRDS(file = here::here("output", "catch_at_length_1530.rds")) |>
  dplyr::ungroup() |>
  dplyr::mutate(CPUE_N_KM2 = FREQ_EXPANDED/AREA_SWEPT_KM2) |>
  dplyr::select(CPUE_N_KM2, TREATMENT, SPECIES_CODE, MATCHUP, SIZE_BIN) |>
  tidyr::pivot_wider(id_cols = c("MATCHUP", "SPECIES_CODE", "SIZE_BIN"), 
                     names_from = "TREATMENT",
                     values_from = "CPUE_N_KM2",
                     values_fill = 0) |>
  dplyr::mutate(obs_ratio = `30`/`15`,
                method = "SCCAL") |>
  dplyr::inner_join(sratio::data_1530$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(), 
                    by = "MATCHUP")

pratio_samples <- readRDS(here::here("output", "pratio_samples.rds")) |>
  dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)),
                method = "SR") |>
  dplyr::inner_join(sratio::data_1530$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(),
                    by = "MATCHUP")

for(ii in 1:length(sp_codes)) {
  
  bootstrap_fit <- data.frame()
  
  sccal_path <- sccal_bootstrap_path[grep(pattern = sp_codes[ii], x = sccal_bootstrap_path)]
  
  sratio_path <- sratio_bootstrap_path[grep(pattern = sp_codes[ii], x = sratio_bootstrap_path)]
  
  if(length(sccal_path) > 0) {
    bootstrap_fit <- dplyr::bind_rows(bootstrap_fit,
                                readRDS(sccal_path) |>
                                  dplyr::mutate(method = "SCCAL"))
  }

  if(length(sratio_path) > 0) {
    bootstrap_fit <- dplyr::bind_rows(bootstrap_fit,
                                readRDS(sratio_path) |>
                                  dplyr::mutate(method = "SR"))
  }
  
  sp_observations <- dplyr::filter(obs_ratio, SPECIES_CODE == sp_codes[ii])
  
  bootstrap_quantiles <- bootstrap_fit |>
    dplyr::group_by(SIZE_BIN, method) |>
    dplyr::summarise(p_q025 = quantile(p12, 0.025),
                     p_q250 = quantile(p12, 0.25),
                     p_q500 = quantile(p12, 0.5),
                     p_q750 = quantile(p12, 0.75),
                     p_q975 = quantile(p12, 0.975),
                     sratio_q025 = quantile(s21, 0.025),
                     sratio_q250 = quantile(s21, 0.25),
                     sratio_q500 = quantile(s21, 0.5),
                     sratio_q750 = quantile(s21, 0.75),
                     sratio_q975 = quantile(s21, 0.975))
  
  plot_obs_histogram <- ggplot() +
    geom_histogram(data = sp_observations,
                   mapping = aes(x = SIZE_BIN, fill = factor(YEAR)),
                   bins = length(unique(sp_observations$SIZE_BIN))-1) +
    scale_x_continuous(name = sratio:::species_code_label(sp_codes[ii]), expand = c(0,0)) +
    scale_y_continuous(name = "Matchups (#)") +
    scale_fill_colorblind() +
    theme_bw() +
    theme(legend.position = c(0.17,0.87),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size = 6.5),
          legend.key.height = unit(2, units = "mm"),
          legend.key.width = unit(4, units = "mm"))
  
  plot_ratios <- ggplot() +
    geom_point(data = sp_observations,
               mapping = aes(x = SIZE_BIN,
                             y = obs_ratio,
                             color = method),
               size = rel(0.3),
               alpha = 0.5) +
    geom_point(data = dplyr::filter(pratio_samples, SPECIES_CODE == sp_codes[ii]),
               mapping = aes(x = SIZE_BIN,
                             y = 1/p-1,
                             color = method),
               size = rel(0.3),
               alpha = 0.5) +
    geom_ribbon(data = bootstrap_quantiles,
                mapping = aes(x = SIZE_BIN, 
                              ymin = sratio_q025, 
                              ymax = sratio_q975,
                              fill = method),
                alpha = 0.2) +
    geom_ribbon(data = bootstrap_quantiles,
                mapping = aes(x = SIZE_BIN, 
                              ymin = sratio_q250, 
                              ymax = sratio_q750,
                              fill = method),
                alpha = 0.4) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN, 
                            y = sratio_q500,
                            color = method)) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_x_continuous(sratio:::species_code_label(x = sp_codes[ii]))  +
    scale_fill_manual(values = c("SCCAL" = "#01665E", "SR" = "#8C510A")) +
    scale_color_manual(values = c("SCCAL" = "#01665E", "SR" = "#8C510A")) +
    scale_y_continuous(name = expression(italic(S['L,15,30'])),
                       limits = c(0, ifelse(max(bootstrap_quantiles$sratio_q975) > 10, 10, 
                                            max(bootstrap_quantiles$sratio_q975)))) +
    theme_bw() +
    theme(legend.position = c(0.83, 0.87), 
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size = 6.5),
          legend.key.height = unit(2, units = "mm"),
          legend.key.width = unit(4, units = "mm"))
  
  png(filename = here::here("plots", paste0(sp_codes[ii], "_sccal_vs_sratio.png")), width = 113, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_ratios,
                           nrow = 1,
                           labels = LETTERS[1:2]))
  dev.off()
  
  sccal_vs_sr <- dplyr::inner_join(
    sp_observations |>
      dplyr::select(SPECIES_CODE, SIZE_BIN, MATCHUP, obs_ratio_sccal = obs_ratio),
    dplyr::filter(pratio_samples, SPECIES_CODE == sp_codes[ii]) |>
      dplyr::mutate(sratio = 1/p-1) |>
      dplyr::select(SPECIES_CODE, SIZE_BIN, MATCHUP, obs_ratio_sr = sratio)
  )
  
  png(filename = here::here("plots", paste0(sp_codes[ii], "_sccal_vs_sratio_observations.png")), width = 70, height = 70, units = "mm", res = 300)
  print(ggplot() +
          geom_abline(intercept = 0, slope = 1) +
          geom_point(data = sccal_vs_sr, mapping = aes(x = obs_ratio_sccal, 
                                                y = 1/obs_ratio_sr)) +
          scale_x_continuous(name = expression(italic(S['L,15,30']~(SCCAL)))) +
          scale_y_continuous(name = expression(italic(S['L,15,30']~(SR)))))
  dev.off()
  
  
}
