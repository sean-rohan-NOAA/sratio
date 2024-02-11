# Plot catch-at-length model results
library(sratio)

bootstrap_results_path <- list.files(here::here("output"), 
                                     recursive = TRUE, 
                                     pattern = "sccal_model_bootstrap_results_", 
                                     full.names = TRUE)

obs_ratio <- readRDS(file = here::here("output", "catch_at_length_1530.rds")) |>
  dplyr::ungroup() |>
  dplyr::mutate(CPUE_N_KM2 = FREQ_EXPANDED/AREA_SWEPT_KM2) |>
  dplyr::select(CPUE_N_KM2, TREATMENT, SPECIES_CODE, MATCHUP, SIZE_BIN) |>
  tidyr::pivot_wider(id_cols = c("MATCHUP", "SPECIES_CODE", "SIZE_BIN"), 
                     names_from = "TREATMENT",
                     values_from = "CPUE_N_KM2",
                     values_fill = 0) |>
  dplyr::mutate(obs_ratio = `30`/`15`)

for(ii in 1:length(bootstrap_results_path)) {
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_results_path[ii])))
  
  boot_fit <- readRDS(file = bootstrap_results_path[ii])
  
  boot_fit_quantile <- dplyr::group_by(boot_fit, SIZE_BIN, TREATMENT) |>
    dplyr::summarise(mean = mean(fit),
                     q025 = quantile(fit, 0.025),
                     q250 = quantile(fit, 0.25),
                     q750 = quantile(fit, 0.75),
                     q975 = quantile(fit, 0.975),
                     .groups = "keep")
  
  y_lim <- quantile(boot_fit$fit, c(0.01, 0.99))
  
  ggplot() +
    geom_ribbon(data = dplyr::filter(boot_fit_quantile, TREATMENT == 15), 
                mapping = aes(x = SIZE_BIN, ymin = q025, ymax = q975, fill = TREATMENT), 
                alpha = 0.2) +
    geom_ribbon(data = dplyr::filter(boot_fit_quantile, TREATMENT == 15), 
                mapping = aes(x = SIZE_BIN, ymin = q250, ymax = q750, fill = TREATMENT), 
                alpha = 0.4) +
    geom_path(data = dplyr::filter(boot_fit_quantile, TREATMENT == 15), 
              mapping = aes(x = SIZE_BIN, y = mean, color = TREATMENT)) +
    geom_ribbon(data = dplyr::filter(boot_fit_quantile, TREATMENT == 30), 
                mapping = aes(x = SIZE_BIN, ymin = q025, ymax = q975, fill = TREATMENT), 
                alpha = 0.2) +
    geom_ribbon(data = dplyr::filter(boot_fit_quantile, TREATMENT == 30), 
                mapping = aes(x = SIZE_BIN, ymin = q250, ymax = q750, fill = TREATMENT), 
                alpha = 0.4) +
    geom_path(data = dplyr::filter(boot_fit_quantile, TREATMENT == 30), 
              mapping = aes(x = SIZE_BIN, y = mean, color = TREATMENT)) +
    scale_x_continuous(sratio:::species_code_label(x = sp_code))  +
    scale_y_continuous(name = expression(CPUE~'#'/km^2)) +
    scale_fill_colorblind() +
    scale_color_colorblind() +
    theme_bw()
  
  
  boot_wide <- tidyr::pivot_wider(boot_fit, 
                                  names_from = TREATMENT, 
                                  values_from = fit, 
                                  values_fill = 0) |>
    as.data.frame()
  
  boot_wide$ratio <- boot_wide[, ncol(boot_wide)]/boot_wide[, (ncol(boot_wide)-1)]
  
  boot_fit_ratio <- boot_wide |>
    dplyr::group_by(SIZE_BIN) |>
    dplyr::summarise(q025 = quantile(ratio, 0.025),
                     q250 = quantile(ratio, 0.25),
                     q500 = quantile(ratio, 0.5),
                     q750 = quantile(ratio, 0.75),
                     q975 = quantile(ratio, 0.975),
                     .groups = "keep")
  
  
  png(filename = here::here("plots", paste0(sp_code, "_total_catch_gam_ratio_ribbon.png")), width = 120, height = 120, units = "mm", res = 300)
  print(
    ggplot() +
      geom_point(data = dplyr::filter(obs_ratio, SPECIES_CODE == sp_code),
                 mapping = aes(x = SIZE_BIN, y = obs_ratio),
                 size = rel(0.3),
                 alpha = 0.5) +
      geom_ribbon(data = boot_fit_ratio,
                  mapping = aes(x = SIZE_BIN, ymin = q025, ymax = q975),
                  alpha = 0.2) +
      geom_ribbon(data = boot_fit_ratio,
                  mapping = aes(x = SIZE_BIN, ymin = q250, ymax = q750),
                  alpha = 0.4) +
      geom_path(data = boot_fit_ratio,
                mapping = aes(x = SIZE_BIN, y = q500)) +
      geom_hline(yintercept = 1, linetype = 2) +
      scale_x_continuous(sratio:::species_code_label(x = sp_code))  +
      scale_y_continuous(name = expression(CPUE['L,30']/CPUE['L,15']),
                         limits = c(0, ifelse(max(boot_fit_ratio$q975) > 10, 10, max(boot_fit_ratio$q975)))) +
      theme_bw()
  )
  dev.off()
  
  # png(filename = here::here("plots", paste0(sp_code, "_total_catch_gam_ratio_lines.png")), width = 120, height = 120, units = "mm", res = 300)
  # print(
  #   boot_fit_ratio |>
  #     tidyr::pivot_longer(cols = c("q025", "q250", "q500", "q750", "q975")) |>
  #     dplyr::inner_join(data.frame(name = c("q025", "q250", "q500", "q750", "q975"),
  #                                  level = c("95%", "50%", "Median", "50%", "95%"))) |>
  #     ggplot() +
  #     geom_hline(yintercept = 1, col = "darkblue") +
  #     geom_path(mapping = aes(x = SIZE_BIN, y = value, linetype = level, group = name)) +
  #     scale_linetype_manual(values = c("95%" = 3, "50%" = 2, "Median" = 1)) +
  #     scale_x_continuous(sratio:::species_code_label(x = sp_code))  +
  #     scale_y_continuous(name = expression(CPUE['L,30']/CPUE['L,15']),
  #                        limits = c(0, ifelse(max(boot_fit_ratio$q975) > 10, 10, max(boot_fit_ratio$q975)))) +
  #     theme_bw() +
  #     theme(legend.position = "none")
  # )
  # dev.off()
  
}
