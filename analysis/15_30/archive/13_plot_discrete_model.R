library(sratio)

bootstrap_results_path <- list.files(here::here("analysis", "15_30", "output"), 
                                     recursive = TRUE, 
                                     pattern = "discrete_bootstrap_results_", 
                                     full.names = TRUE)

pratio_samples <- readRDS(here::here("analysis", "15_30", "output", "discrete_samples.rds")) |>
  dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP)),
                SIZE_BIN = as.numeric(as.character(SIZE_BIN))) |>
  dplyr::inner_join(sratio::data_1530$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(),
                    by = "MATCHUP")

dir.create(here::here("analysis", "15_30",  
                      "plots", "discrete_fit"),
           showWarnings = FALSE)

for(ii in 1:length(bootstrap_results_path)) {
  
  bootstrap_df <- readRDS(file = bootstrap_results_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_results_path[ii])))
  
  bootstrap_quantiles <- bootstrap_df |>
    dplyr::group_by(SIZE_BIN) |>
    dplyr::summarise(p_q025 = quantile(p12, 0.025, na.rm = TRUE),
                     p_q250 = quantile(p12, 0.25, na.rm = TRUE),
                     p_q500 = quantile(p12, 0.5, na.rm = TRUE),
                     p_q750 = quantile(p12, 0.75, na.rm = TRUE),
                     p_q975 = quantile(p12, 0.975, na.rm = TRUE),
                     sratio_q025 = quantile(s12, 0.025, na.rm = TRUE),
                     sratio_q250 = quantile(s12, 0.25, na.rm = TRUE),
                     sratio_q500 = quantile(s12, 0.5, na.rm = TRUE),
                     sratio_q750 = quantile(s12, 0.75, na.rm = TRUE),
                     sratio_q975 = quantile(s12, 0.975, na.rm = TRUE)) |>
    dplyr::mutate(type = "Bootstrap")
  
  
  # Make plots of catch ratio and selectivity ratio ----
  plot_pratio <- ggplot() +
    geom_point(data = dplyr::filter(pratio_samples, 
                                    SPECIES_CODE == sp_code),
               mapping = aes(x = SIZE_BIN, y = p),
               size = rel(0.3),
               alpha = 0.5) +
    geom_ribbon(data = bootstrap_quantiles,
                mapping = aes(x = SIZE_BIN,
                              ymin = p_q025,
                              max = p_q975),
                alpha = 0.5,
                fill = "grey20") +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            y = p_q250),
              linetype = 3) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            y = p_q750),
              linetype = 3) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            y = p_q500)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    scale_x_continuous(name = sratio:::species_code_label(sp_code)) +
    scale_y_continuous(name = expression(italic(p['L,30,15'])), limits = c(0,1)) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  hist_df <- dplyr::filter(pratio_samples, SPECIES_CODE == sp_code)
  
  plot_obs_histogram <- ggplot() +
    geom_histogram(data = hist_df ,
                   mapping = aes(x = SIZE_BIN, fill = factor(YEAR)),
                   bins = length(unique(hist_df $SIZE_BIN))-1) +
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
  
  plot_sratio <- ggplot() +
    geom_point(data = dplyr::filter(pratio_samples, SPECIES_CODE == sp_code),
               mapping = aes(x = SIZE_BIN, y = p/(1-p)),
               size = rel(0.3),
               alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_ribbon(data = bootstrap_quantiles,
                mapping = aes(x = SIZE_BIN,
                              ymin = sratio_q025,
                              max = sratio_q975),
                alpha = 0.5,
                fill = "grey20") +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            y = sratio_q250),
              linetype = 3) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            y = sratio_q750),
              linetype = 3) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            y = sratio_q500)) +
    scale_x_continuous(name = sratio:::species_code_label(sp_code)) +
    scale_y_log10(name = expression(italic(S['L,30,15'])~(SR)), expand = c(0.05, 0.05)) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  # Write plots to file
  ragg::agg_png(file = here::here("analysis", "15_30", 
                                  "plots", paste0(sp_code, "_discrete_three_panel.png")), 
                width = 169, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_pratio,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:3]))
  dev.off()
  
  ragg::agg_png(file = here::here("analysis", "15_30",  
                                  "plots", "discrete_fit",
                                  paste0(sp_code, "_discrete_two_panel.png")), 
                width = 104, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:3]))
  dev.off()
  
}
