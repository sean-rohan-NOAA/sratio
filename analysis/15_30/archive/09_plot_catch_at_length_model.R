# Plot catch-at-length model results
library(sratio)
library(shadowtext)

bootstrap_results_path <- list.files(here::here("analysis", "15_30", "output"), 
                                     recursive = TRUE, 
                                     pattern = "sccal_bootstrap_results_", 
                                     full.names = TRUE)

dir.create(here::here("analysis", "15_30",  
                      "plots", "sccal_fit"),
           showWarnings = FALSE)

obs_ratio <- readRDS(file = here::here("analysis", "15_30", "output", "catch_at_length_1530.rds")) |>
  dplyr::ungroup() |>
  dplyr::mutate(CPUE_N_KM2 = FREQ_EXPANDED/AREA_SWEPT_KM2) |>
  dplyr::select(CPUE_N_KM2, TREATMENT, SPECIES_CODE, MATCHUP, SIZE_BIN) |>
  tidyr::pivot_wider(id_cols = c("MATCHUP", "SPECIES_CODE", "SIZE_BIN"), 
                     names_from = "TREATMENT",
                     values_from = "CPUE_N_KM2",
                     values_fill = 0) |>
  dplyr::mutate(obs_ratio = `30`/`15`) |>
  dplyr::inner_join(sratio::data_1530$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(), 
                    by = "MATCHUP") |>
  dplyr::mutate(COMMON_NAME = sratio::species_code_label(x = SPECIES_CODE, 
                                                         type = "common_name", 
                                                         make_factor = TRUE))

year_colors <- c(`1995` = "#0072B2", 
                 `1998` =  "#F0E442", 
                 `2021` =  "#009E73", 
                 `2022` =  "#56B4E9", 
                 `2023` = "#000000", 
                 `2024` = "#E69F00")

fits_all_species <- data.frame()

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
  
  fits_all_species <- dplyr::bind_rows(fits_all_species, boot_fit_quantile)
  
  y_lim <- quantile(boot_fit$s12, c(0.01, 0.99))
  
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
    scale_y_continuous(name = expression(italic(S['L,30,15']~(SCCAL))),
                       limits = c(0, ifelse(max(boot_fit_quantile$sratio_q975) > 10, 10, 
                                            max(boot_fit_quantile$sratio_q975)))) +
    theme_bw()
  
  png(filename = here::here("analysis", "15_30", "plots", "sccal_fit",
                            paste0(sp_code, "_sccal_ratio.png")), 
      width = 113, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_sccal,
                           nrow = 1,
                           labels = LETTERS[1:2]))
  dev.off()
  
}

# Multi-panel plot with all species

fits_all_species <- fits_all_species |>
  dplyr::mutate(COMMON_NAME = sratio::species_code_label(x = SPECIES_CODE, 
                                                         type = "common_name", 
                                                         make_factor = TRUE))

annotate_pollock_30 <- data.frame(x = -Inf, y = Inf, label = "30 higher", 
                                  COMMON_NAME = species_code_label(21740, type = "common_name", make_factor = TRUE))
annotate_pollock_15 <- data.frame(x = -Inf, y = -Inf, label = "15 higher", 
                                  COMMON_NAME = species_code_label(21740, type = "common_name", make_factor = TRUE))

set_y_max <- 4

plot_multipanel_sratio <- ggplot() +
  geom_point(data = obs_ratio,
             mapping = aes(x = SIZE_BIN, y = obs_ratio),
             alpha = 0.2,
             size = rel(0.9),
             color = "grey20") +
  geom_hline(yintercept = 1, linetype = 2,
             linewidth = rel(1),
             color = "red") +
  geom_ribbon(data = 
                dplyr::mutate(
                  fits_all_species, 
                  sratio_q975 = dplyr::if_else(sratio_q975 > set_y_max, set_y_max, sratio_q975)),
              mapping = aes(x = SIZE_BIN,
                            ymin = sratio_q025,
                            max = sratio_q975),
              alpha = 0.5,
              fill = "grey50") +
  geom_path(data = fits_all_species,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q250),
            linetype = 3,
            linewidth = rel(1)) +
  geom_path(data = fits_all_species,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q750),
            linetype = 3,
            linewidth = rel(1)) +
  geom_path(data = fits_all_species,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q500),
            linewidth = rel(1)) +
  geom_shadowtext(data = annotate_pollock_30,
                  mapping = aes(x = x, y = y, label = label), 
                  hjust = -0.1, 
                  vjust = 1.2, 
                  color = "red", 
                  bg.color = "white",
                  fontface = "bold",
                  size = 5) +
  geom_shadowtext(data = annotate_pollock_15,
                  mapping = aes(x = x, y = y, label = label), 
                  hjust = -0.1, 
                  vjust = -0.7, 
                  color = "red", 
                  bg.color = "white",
                  fontface = "bold",
                  size = 5) +
  facet_wrap(~COMMON_NAME, scales = "free", ncol = 3) +
  scale_x_continuous(name = "Size") +
  scale_y_continuous(name = expression(italic(S['L,30,15'])), expand = c(0,0), limits = c(-0.05, set_y_max)) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, units = "mm"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        strip.text = element_text(size = 10))


ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "sccal_fit", "sccal_all_species.png"), 
              width = 169, 
              height = 169*1.25,
              units = "mm",
              res = 300)
print(plot_multipanel_sratio)
dev.off()

