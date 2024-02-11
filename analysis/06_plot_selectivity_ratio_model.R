library(sratio)

bootstrap_results_path <- list.files(here::here("output"), 
                                    recursive = TRUE, 
                                    pattern = "sratio_bootstrap_results_", 
                                    full.names = TRUE)

pratio_samples <- readRDS(here::here("output", "pratio_samples.rds"))

for(ii in 1:length(bootstrap_results_path)) {
  
  bootstrap_df <- readRDS(file = bootstrap_results_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_results_path[ii])))
  
  bootstrap_quantiles <- bootstrap_df |>
    dplyr::group_by(SIZE_BIN) |>
    dplyr::summarise(q025 = quantile(fit, 0.025),
                     q975 = quantile(fit, 0.975),
                     q500 = quantile(fit, 0.500),
                     q250 = quantile(fit, 0.25),
                     q750 = quantile(fit, 0.75)) |> 
    dplyr::mutate(p_q025 = inv_logit(q025),
                  sratio_q025 = 1/inv_logit(q025)-1,
                  p_q975 = inv_logit(q975),
                  sratio_q975 = 1/inv_logit(q975)-1,
                  p_q500 = inv_logit(q500),
                  sratio_q500 = 1/inv_logit(q500)-1,
                  p_q250 = inv_logit(q250),
                  sratio_q250 = 1/inv_logit(q250)-1,
                  p_q750 = inv_logit(q750),
                  sratio_q750 = 1/inv_logit(q750)-1) |>
    dplyr::mutate(type = "Observed")
  
  
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
    geom_histogram(data = hist_df,
                   mapping = aes(x = SIZE_BIN),
                   bins = length(unique(hist_df$SIZE_BIN))-1,
                   fill = "grey50",
                   color = "grey50") +
    scale_x_continuous(name = sratio:::species_code_label(sp_code), expand = c(0,0)) +
    scale_y_continuous(name = "Matchups (#)") +
    theme_bw()
  
  plot_sratio <- ggplot() +
    geom_point(data = dplyr::filter(pratio_samples, SPECIES_CODE == sp_code),
               mapping = aes(x = SIZE_BIN, y = 1/p-1),
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
    # scale_y_continuous(name = expression(italic(S['L,15,30'])), expand = c(0.05, 0.05)) +
    scale_y_log10(name = expression(italic(S['L,15,30'])), expand = c(0.05, 0.05)) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  # Write plots to file
  ragg::agg_png(file = here::here("plots", paste0(sp_code, "_trawl_height_two_panel_ratios_n.png")), width = 113, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_pratio,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:2]))
  dev.off()
  
  # Write plots to file
  ragg::agg_png(file = here::here("plots", paste0(sp_code, "_trawl_height_three_panel_ratios_n.png")), width = 169, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_pratio,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:3]))
  dev.off()
  
}

