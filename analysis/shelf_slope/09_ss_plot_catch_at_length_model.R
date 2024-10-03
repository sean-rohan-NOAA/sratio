# Plot catch-at-length model results
library(sratio)

bootstrap_results_path <- list.files(here::here("analysis", "shelf_slope", "output"), 
                                     recursive = TRUE, 
                                     pattern = "sccal_bootstrap_results_", 
                                     full.names = TRUE)

obs_ratio <- readRDS(file = here::here("analysis", "shelf_slope", "output", "catch_at_length_ss.rds")) |>
  dplyr::ungroup() |>
  dplyr::mutate(CPUE_N_KM2 = FREQ_EXPANDED/AREA_SWEPT_KM2) |>
  dplyr::select(CPUE_N_KM2, TREATMENT, SPECIES_CODE, MATCHUP, SIZE_BIN) |>
  tidyr::pivot_wider(id_cols = c("MATCHUP", "SPECIES_CODE", "SIZE_BIN"), 
                     names_from = "TREATMENT",
                     values_from = "CPUE_N_KM2",
                     values_fill = 0) |>
  dplyr::mutate(obs_ratio = `172`/`44`) |>
  dplyr::inner_join(sratio::data_ss$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(), 
                    by = "MATCHUP")

for(ii in 1:length(bootstrap_results_path)) {
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_results_path[ii])))
  
  sp_observations <- dplyr::filter(obs_ratio, SPECIES_CODE == sp_code)
  
  boot_fit <- readRDS(file = bootstrap_results_path[ii])
  
  boot_fit_quantile <- dplyr::group_by(boot_fit, SIZE_BIN, SPECIES_CODE) |>
    dplyr::summarise(sratio_q025 = quantile(s12, 0.025),
                     sratio_q250 = quantile(s12, 0.25),
                     sratio_q500 = quantile(s12, 0.5),
                     sratio_q750 = quantile(s12, 0.75),
                     sratio_q975 = quantile(s12, 0.975),
                     .groups = "keep")
  
  y_lim <- quantile(boot_fit$s21, c(0.01, 0.99))
  
  plot_obs_histogram <- ggplot() +
    geom_histogram(data = sp_observations,
                   mapping = aes(x = SIZE_BIN, fill = factor(YEAR)),
                   bins = length(unique(sp_observations$SIZE_BIN))-1) +
    scale_x_continuous(name = sratio:::species_code_label(sp_code), expand = c(0,0)) +
    scale_y_continuous(name = "Matchups (#)") +
    scale_fill_colorblind() +
    theme_bw() +
    theme(legend.position = c(0.17,0.87),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size = 6.5),
          legend.key.height = unit(2, units = "mm"),
          legend.key.width = unit(4, units = "mm"))
  
  plot_sccal <-  ggplot() +
    geom_point(data = sp_observations,
               mapping = aes(x = SIZE_BIN, 
                             y = obs_ratio),
               size = rel(0.3),
               alpha = 0.5) +
    geom_ribbon(data = boot_fit_quantile,
                mapping = aes(x = SIZE_BIN, 
                              ymin = sratio_q025, 
                              ymax = sratio_q975),
                alpha = 0.2) +
    geom_ribbon(data = boot_fit_quantile,
                mapping = aes(x = SIZE_BIN, 
                              ymin = sratio_q250, 
                              ymax = sratio_q750),
                alpha = 0.4) +
    geom_path(data = boot_fit_quantile,
              mapping = aes(x = SIZE_BIN, 
                            y = sratio_q500)) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_x_continuous(sratio:::species_code_label(x = sp_code))  +
    scale_y_continuous(name = expression(italic(S['L,PNE,83-112']~(SCCAL))),
                       limits = c(0, ifelse(max(boot_fit_quantile$sratio_q975) > 10, 10, 
                                            max(boot_fit_quantile$sratio_q975)))) +
    theme_bw()
  
  png(filename = here::here("analysis", "shelf_slope", "plots", 
                            paste0(sp_code, "_sccal_ratio.png")), 
      width = 113, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_sccal,
                           nrow = 1,
                           labels = LETTERS[1:2]))
  dev.off()
  
}
