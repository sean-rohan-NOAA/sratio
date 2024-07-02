library(akgfmaps)
library(ggthemes)

map_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")

ggplot() +
  geom_sf(data = map_layers$survey.grid) +
  geom_sf_text(data = sf::st_centroid(map_layers$survey.grid),
          mapping = aes(label = STATIONID),
          size = 2.2) +
  geom_sf(data = map_layers$akland) +
  scale_x_continuous(limits = map_layers$plot.boundary$x) +
  scale_y_continuous(limits = map_layers$plot.boundary$y) +
  theme_bw()

area_4 <- dplyr::filter(map_layers$survey.grid, 
                        STATIONID %in% c(paste0("S-", 22:28),
                                         paste0("T-", 25:28),
                                         paste0("U-", 25:28),
                                         paste0("V-", 25:28))
) |>
  dplyr::mutate(area = "Area 4 (20%)",
                area_n = "Area 4 (4 stations)")

area_3 <- dplyr::filter(map_layers$survey.grid, 
                        STATIONID %in% c(paste0("Q-", 26:31),
                                         paste0("P-", 26:32),
                                         paste0("O-", 26:31),
                                         paste0("N-", 26:31), 
                                         "QP2726", "PO2726")
) |>
  dplyr::mutate(area = "Area 3 (20%)",
                area_n = "Area 3 (4 stations)")

area_2 <- dplyr::filter(map_layers$survey.grid, 
                        STATIONID %in% c("E-22", 
                                         paste0("F-", 22:25),
                                         paste0("G-", 24:26),
                                         "H-25", "H-26", "I-25", "I-26", "J-25", "J-26")
) |>
  dplyr::mutate(area = "Area 2 (20%)",
                area_n = "Area 2 (4 stations)")

area_1 <- dplyr::filter(map_layers$survey.grid, 
                        STATIONID %in% c(paste0("K-", 10:14),
                                         paste0("J-", 10:14),
                                         paste0("I-", 10:14),
                                         paste0("H-", 10:14),
                                         paste0("G-", 10:14),
                                         paste0("F-", 10:14),
                                         paste0("E-", 10:12),
                                         paste0("D-", 10),
                                         paste0("K-0", 8:9),
                                         paste0("J-0", 8:9),
                                         paste0("I-0", 8:9),
                                         paste0("H-0", 8:9),
                                         paste0("G-0", 8:9),
                                         paste0("F-0", 8:9),
                                         paste0("E-0", 8:9),
                                         paste0("D-0", 8:9))
) |>
  dplyr::mutate(area = "Area 1 (40%)",
                area_n = "Area 1 (8 stations)")

sample_zones <- dplyr::bind_rows(area_1, area_2, area_3, area_4) |>
  dplyr::group_by(area, area_n) |>
  summarise(do_union = TRUE)

sample_zones |>
  sf::st_transform(crs = "WGS84") |>
  sf::st_make_valid() |>
  sf::st_write(here::here("analysis", "15_30", "output", "2024_samples_zones_15_30.shp"),
               append = FALSE)

ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "map_target_zones_2024.png"), 
              width = 6, height = 4, units = "in", res = 300)
print(
ggplot() +
  geom_sf(data = map_layers$survey.grid, fill = NA) +
  geom_sf(data = sample_zones, 
          mapping = aes(fill = area),
          alpha = 0.7) +
  geom_sf(data = map_layers$akland, 
          linewidth = 0.2, 
          fill = "grey50", 
          color = "black") +
  scale_x_continuous(limits = map_layers$plot.boundary$x, 
                     breaks = map_layers$lon.breaks) +
  scale_y_continuous(limits = map_layers$plot.boundary$y, 
                     breaks = map_layers$lat.breaks) +
  scale_fill_colorblind(name = "Area") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.17),
        legend.title = element_blank(),
        legend.text = element_text(size = 6.8))
)
dev.off()


ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "map_target_zones_n_2024.png"), 
              width = 6, height = 4, units = "in", res = 300)
print(
  ggplot() +
    geom_sf(data = map_layers$survey.grid, fill = NA) +
    geom_sf(data = sample_zones, 
            mapping = aes(fill = area_n),
            alpha = 0.7) +
    geom_sf(data = map_layers$akland, 
            linewidth = 0.2, 
            fill = "grey50", 
            color = "black") +
    scale_x_continuous(limits = map_layers$plot.boundary$x, 
                       breaks = map_layers$lon.breaks) +
    scale_y_continuous(limits = map_layers$plot.boundary$y, 
                       breaks = map_layers$lat.breaks) +
    scale_fill_colorblind(name = "Area") +
    theme_bw() +
    theme(legend.position = c(0.15, 0.17),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.background = element_blank())
)
dev.off()

