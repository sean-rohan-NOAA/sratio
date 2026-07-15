library(akgfmaps) # github..com/afsc-gap-products/akgfmaps
library(shadowtext)
library(ggthemes)

ebs_layers <- get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")

haul_locs <- 
  sratio::data_1530$haul |>
  dplyr::filter(TREATMENT == 15) |>
  sf::st_as_sf(
    coords = c("START_LONGITUDE", "START_LATITUDE"),
    crs = "WGS84"
  ) |>
  sf::st_transform(crs = "EPSG:3338") |>
  dplyr::arrange(-YEAR)

bathymetry_labels <- 
  data.frame(depth = c("50 m", "100 m", "200 m"),
             x = c(-165, -168, -170),
             y = c(57.8, 56.7, 56)) |>
  sf::st_as_sf(
    coords = c("x", "y"),
    crs = "WGS84"
  ) |>
  sf::st_transform(crs = "EPSG:3338")

bathymetry_labels[c("x", "y")] <- sf::st_coordinates(bathymetry_labels)


land_label <-
  data.frame(
    label = "Alaska",
    x = -160,
    y = 62
  ) |>
  sf::st_as_sf(
    coords = c("x", "y"),
    crs = "WGS84"
  ) |>
  sf::st_transform(crs = "EPSG:3338")

land_label[c("x", "y")] <- sf::st_coordinates(land_label)


water_labels <-
  data.frame(
    label = "Bristol Bay",
    x = -161,
    y = 57.5
  ) |>
  sf::st_as_sf(
    coords = c("x", "y"),
    crs = "WGS84"
  ) |>
  sf::st_transform(crs = "EPSG:3338")

water_labels[c("x", "y")] <- sf::st_coordinates(water_labels)

# 2. Create the map using ggplot2
p_map <- 
  ggplot() +
  geom_sf(data = ebs_layers$survey.area, fill = NA, color = "black", alpha = 0.7) +
  geom_sf(data = ebs_layers$akland, fill = "grey80", color = "grey40") +
  geom_sf(data = ebs_layers$bathymetry, color = "skyblue", alpha = 0.5, linewidth = 0.3) +
  geom_sf(
    data = haul_locs,
          mapping = aes(color = factor(YEAR),
                        shape = factor(YEAR)),
          size = 1.5
  ) +
  geom_shadowtext(
    data = bathymetry_labels,
    mapping = aes(x = x, y = y, label = depth),
    bg.color = "white",
    color = "skyblue",
    size = 3
  ) +
  geom_shadowtext(
    data = land_label,
    mapping = aes(x = x, y = y, label = label),
    bg.color = "white",
    color = "black",
    size = 4
  ) +
  # geom_shadowtext(
  #   data = water_labels,
  #   mapping = aes(x = x, y = y, label = label),
  #   bg.color = "white",
  #   color = "black",
  #   size = 2.5
  # ) +
  scale_x_continuous(
    limits = ebs_layers$plot.boundary$x,
    breaks = ebs_layers$lon.breaks
  ) +
  scale_y_continuous(
    limits = ebs_layers$plot.boundary$y,
    breaks = ebs_layers$lat.breaks
  ) +
  scale_shape(name = "Year", solid = FALSE) +
  scale_color_tableau(name = "Year") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.12,0.21),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.title = element_blank(),
    legend.key.height = unit(3, units = "mm"),
    legend.key.width = unit(3, units = "mm"),
    panel.grid = element_line(color = "grey80", linewidth = 0.2)
  )

png(filename = here::here("analysis", "somerton_2002", "plots", "map_15_30_hauls.png"), width = 100, height = 80, units = "mm", res = 600)
print(p_map)
dev.off()

