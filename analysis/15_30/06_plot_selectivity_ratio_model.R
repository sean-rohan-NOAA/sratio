library(sratio)

bootstrap_results_path <- list.files(here::here("analysis", "15_30", "output"), 
                                    recursive = TRUE, 
                                    pattern = "sratio_bootstrap_results_", 
                                    full.names = TRUE)

pratio_samples <- readRDS(here::here("analysis", "15_30", "output", "pratio_samples.rds")) |>
  dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP))) |>
  dplyr::inner_join(sratio::data_1530$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(),
                    by = "MATCHUP") |>
  dplyr::mutate(COMMON_NAME = sratio::species_code_label(x = SPECIES_CODE, 
                                                         type = "common_name", 
                                                         make_factor = TRUE))

dir.create(here::here("analysis", "15_30",  
                      "plots", "sratio_fit"),
           showWarnings = FALSE)

year_colors <- c(`1995` = "#0072B2", 
                 `1998` =  "#F0E442", 
                 `2021` =  "#009E73", 
                 `2022` =  "#56B4E9", 
                 `2023` = "#000000", 
                 `2024` = "#E69F00")

# set_species_order <- function(SPECIES_CODE, type = "common_name") {
#   
#   output <- factor(SPECIES_CODE, 
#          levels = c(21740,
#                     21720,
#                     11111,
#                     10115,
#                     10110,
#                     10115,
#                     10210,
#                     10261,
#                     10130,
#                     471,
#                     69322,
#                     68580,
#                     68560),
#          labels = c("walleye pollock", 
#                     "Pacific cod",
#                     "sablefish",
#                     "Greenland turbot",
#                     "arrowtooth flounder", 
#                     "Kamchatka flounder",
#                     "yellowfin sole",
#                     "northern rock sole", 
#                     "flathead sole", 
#                     "Alaska skate",
#                     "red king crab",
#                     "snow crab", 
#                     "Tanner crab"))
# 
#   return(output)
#   
# }

fits_all_species <- data.frame()

for(ii in 1:length(bootstrap_results_path)) {
  
  bootstrap_df <- readRDS(file = bootstrap_results_path[ii])
  
  sp_code <- as.numeric(gsub("[^0-9]", "", basename(bootstrap_results_path[ii])))
  
  bootstrap_quantiles <- bootstrap_df |>
    dplyr::group_by(SIZE_BIN, SPECIES_CODE) |>
    dplyr::summarise(p_q025 = quantile(p12, 0.025),
                     p_q250 = quantile(p12, 0.25),
                     p_q500 = quantile(p12, 0.5),
                     p_q750 = quantile(p12, 0.75),
                     p_q975 = quantile(p12, 0.975),
                     sratio_q025 = quantile(s12, 0.025),
                     sratio_q250 = quantile(s12, 0.25),
                     sratio_q500 = quantile(s12, 0.5),
                     sratio_q750 = quantile(s12, 0.75),
                     sratio_q975 = quantile(s12, 0.975)) |>
    dplyr::mutate(type = "Bootstrap")
  
  fits_all_species <- dplyr::bind_rows(fits_all_species, bootstrap_quantiles)
  
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
    geom_histogram(data = dplyr::select(hist_df, MATCHUP, SIZE_BIN, YEAR) |>
                     unique(),
                   mapping = aes(x = SIZE_BIN, fill = factor(YEAR)),
                   bins = length(unique(hist_df $SIZE_BIN))-1) +
    scale_x_continuous(name = sratio:::species_code_label(sp_code), expand = c(0,0)) +
    scale_y_continuous(name = "Matchups (#)") +
    scale_fill_manual(values = year_colors) +
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
                                  "plots", "sratio_fit", paste0(sp_code, "_sratio_three_panel.png")), 
                width = 169, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_pratio,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:3]))
  dev.off()
  
  ragg::agg_png(file = here::here("analysis", "15_30",  
                                  "plots", "sratio_fit",
                                  paste0(sp_code, "_sratio_two_panel.png")), 
                width = 104, height = 70, units = "mm", res = 300)
  print(cowplot::plot_grid(plot_obs_histogram,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:3]))
  dev.off()
  
}


# Multi-panel plot with all species

fits_all_species <- fits_all_species |>
  dplyr::mutate(COMMON_NAME = sratio::species_code_label(x = SPECIES_CODE, 
                                                         type = "common_name", 
                                                         make_factor = TRUE))

plot_multipanel_sratio <- ggplot() +
  geom_point(data = pratio_samples,
             mapping = aes(x = SIZE_BIN, y = p/(1-p)),
             alpha = 0.2,
             size = rel(0.9)) +
  geom_hline(yintercept = 1, linetype = 2,
             linewidth = rel(1),
             color = "red") +
  geom_ribbon(data = fits_all_species,
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
  facet_wrap(~COMMON_NAME, scales = "free", ncol = 3) +
  scale_x_continuous(name = "Size") +
  scale_y_continuous(name = expression(italic(S['L,30,15'])), expand = c(0,0), limits = c(-0.05, 3)) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, units = "mm"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text = element_text(size = 10))


ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "sratio_all_species.png"), 
              width = 169, 
              height = 169*1.25,
              units = "mm",
              res = 300)
print(plot_multipanel_sratio)
dev.off()
