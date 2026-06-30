# Trawl gear geometry plots and summary statistics

library(sratio)

haul_geometry <- 
  sratio::data_ss$haul |>
  dplyr::mutate(WIDTH_TO_HEIGHT = NET_WIDTH/NET_HEIGHT) |>
  dplyr::select(TREATMENT, MATCHUP, NET_WIDTH, NET_HEIGHT, WIDTH_TO_HEIGHT) |>
  tidyr::pivot_wider(id_cols = "MATCHUP", names_from = "TREATMENT", values_from = c("NET_WIDTH", "NET_HEIGHT", "WIDTH_TO_HEIGHT"))

mean(haul_geometry$NET_WIDTH_44 - haul_geometry$NET_WIDTH_172, na.rm = TRUE)
mean(haul_geometry$NET_HEIGHT_44 - haul_geometry$NET_HEIGHT_172, na.rm = TRUE)


net_mean <- 
  dplyr::group_by(sratio::data_ss$haul, TREATMENT) |>
  dplyr::summarise(
    MEAN_NET_WIDTH = mean(NET_WIDTH, na.rm = TRUE),
    MIN_NET_WIDTH = min(NET_WIDTH, na.rm = TRUE),
    MAX_NET_WIDTH = max(NET_WIDTH, na.rm = TRUE),
    MEAN_NET_HEIGHT = mean(NET_HEIGHT, na.rm = TRUE),
    MIN_NET_HEIGHT = min(NET_HEIGHT, na.rm = TRUE),
    MAX_NET_HEIGHT = max(NET_HEIGHT, na.rm = TRUE),
    MEAN_WIDTH_TO_HEIGHT = mean(NET_WIDTH/NET_HEIGHT, na.rm = TRUE)
  )

p_gear_geometry <-
  cowplot::plot_grid(
  ggplot() +
    geom_point(
      data = haul_geometry,
      mapping = aes(x = NET_WIDTH_44, y = NET_WIDTH_172)
    ) +
    geom_hline(
      data = net_mean[net_mean$TREATMENT == "172", ],
      mapping = aes(yintercept = MEAN_NET_WIDTH),
      linetype = 2
               ) +
    geom_vline(
      data = net_mean[net_mean$TREATMENT == "44", ],
      mapping = aes(xintercept = MEAN_NET_WIDTH),
      linetype = 2
    ) +
    scale_x_continuous(name = "83-112 spread (m)") +
    scale_y_continuous(name = "PNE spread (m)") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 8)),
  ggplot() +
    geom_point(
      data = haul_geometry,
      mapping = aes(x = NET_HEIGHT_44, y = NET_HEIGHT_172)
    ) +
    geom_hline(
      data = net_mean[net_mean$TREATMENT == "172", ],
      mapping = aes(yintercept = MEAN_NET_HEIGHT),
      linetype = 2
    ) +
    geom_vline(
      data = net_mean[net_mean$TREATMENT == "44", ],
      mapping = aes(xintercept = MEAN_NET_HEIGHT),
      linetype = 2
    ) +
    scale_x_continuous(name = "83-112 height (m)") +
    scale_y_continuous(name = "PNE height (m)") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 8)),
  ggplot() +
    geom_point(
      data = haul_geometry,
      mapping = aes(x = WIDTH_TO_HEIGHT_44, y = WIDTH_TO_HEIGHT_172)
    ) +
    geom_hline(
      data = net_mean[net_mean$TREATMENT == "172", ],
      mapping = aes(yintercept = MEAN_WIDTH_TO_HEIGHT),
      linetype = 2
    ) +
    geom_vline(
      data = net_mean[net_mean$TREATMENT == "44", ],
      mapping = aes(xintercept = MEAN_WIDTH_TO_HEIGHT),
      linetype = 2
    ) +
    scale_x_continuous(name = "83-112 spread/height") +
    scale_y_continuous(name = "PNE spread/height") +
    theme_bw() +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 8)),
  nrow = 1,
  align = "hv"
)

png(filename = here::here("analysis", "shelf_slope_gear_selectivity", "plots", "gear_geometry.png"),
    width = 169, height = 40, units = "mm", res = 300)
print(p_gear_geometry)
dev.off()



