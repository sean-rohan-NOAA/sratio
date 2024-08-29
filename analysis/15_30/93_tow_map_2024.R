# Locations of 15/30 tows in 2024

library(sratio)
library(akgfmaps)
library(shadowtext)

map_layers <- akgfmaps::get_base_layers(select.region = "sebs", 
                                        set.crs = "EPSG:3338", 
                                        include.corners = FALSE,
                                        use.survey.bathymetry = TRUE)

tows_1530_2024 <- sratio::data_1530$haul |>
  dplyr::filter(CRUISE == 202401) |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"),
                   crs = "EPSG:4269") |>
  dplyr::mutate(label = "15/30 Haul Duration\nComparison Tow")

tow_map_2024 <- ggplot() +
  geom_sf(data = map_layers$akland) +
  geom_sf(data = map_layers$survey.area, 
          fill = NA, 
          color = "black", 
          linewidth = 0.5) +
  geom_sf(data = map_layers$bathymetry,
          color = "grey50") +
  geom_sf(data = tows_1530_2024,
          mapping = aes(color = label),
          size = 3.5) + 
  ggplot2::geom_text(data = subset(map_layers$place.labels, 
                                   type == "mainland"), 
                     aes(x = x, y = y, label = lab), 
                     size = 8, group = 99) + 
  shadowtext::geom_shadowtext(data = subset(map_layers$place.labels, 
                                            type == "peninsula"), 
                              aes(x = x, y = y-1e4, label = lab), size = 5, angle = 30, 
                              bg.color = "white", color = "black", group = 99) + 
  shadowtext::geom_shadowtext(
    data = subset(map_layers$place.labels, 
                  type == "bathymetry"),
    aes(x = x, y = y, label = lab), 
    bg.color = "white", color = "black", 
    size = 4, group = 99) +
  scale_x_continuous(limits = map_layers$plot.boundary$x+c(-5e3, 0),
                     breaks = map_layers$lon.breaks) +
  scale_y_continuous(limits = map_layers$plot.boundary$y,
                     breaks = map_layers$lat.breaks) +
  scale_color_manual(values = "#3EA2D5") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", colour = "grey20"), 
        strip.background = element_rect(fill = "grey85",colour = "grey20"), 
        legend.spacing.y = unit(-0.35, "cm"),
        legend.text = element_text(size = 13),
        legend.background=element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = c(.16, .08),
        legend.box.just = "left",
        legend.box = "vertical", 
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_blank())

ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "sample_map_2024.png"),
              width = 6.5,
              height = 6,
              units = "in",
              res = 300)
print(tow_map_2024)
dev.off()
