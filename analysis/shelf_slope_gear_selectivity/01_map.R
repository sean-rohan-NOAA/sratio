library(akgfmaps) # afsc-gap-products/akgfmaps
library(shadowtext)
library(tidyterra)

station_loc <- 
  sratio::data_ss$haul |>
  dplyr::group_by(YEAR, STATIONID, MATCHUP) |>
  dplyr::summarise(
    START_LONGITUDE = mean(START_LONGITUDE), 
    START_LATITUDE = mean(START_LATITUDE)
  ) |>
  dplyr::mutate(group = ifelse(START_LATITUDE > 57, "North", "South")) |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"), crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338")

bathy_raster <- terra::rast(system.file("./extdata/bathymetry.gpkg", package = "akgfmaps"))

map_layers <- akgfmaps::get_base_layers(select.region = c("sebs", "ebs.slope"), set.crs = "EPSG:3338")

map_layers$survey.area <- map_layers$survey.area |> 
  dplyr::mutate(
  SURVEY_ABBV = ifelse(SURVEY_DEFINITION_ID == 98, "EBS Shelf", "EBS Slope") 
    )

area_label <- navmaps::st_primary_point_on_surface(map_layers$survey.area)
area_label[c("x", "y")] <- sf::st_coordinates(area_label)

north_area <- station_loc |>
  dplyr::filter(group == "North") |>
  sf::st_buffer(dist = 2.5e4, endCapStyle = "SQUARE", joinStyle = "MITRE") |>
  sf::st_bbox() |>
  sf::st_as_sfc() |>
  sf::st_as_sf()

bbox_north <- sf::st_bbox(north_area)

south_area <- station_loc |>
  dplyr::filter(group == "South") |>
  sf::st_buffer(dist = 2.5e4, endCapStyle = "SQUARE", joinStyle = "MITRE") |>
  sf::st_bbox() |>
  sf::st_as_sfc() |>
  sf::st_as_sf()

bbox_south <- sf::st_bbox(south_area)

# Create 100 m resolution bathymetry
bathy_south <- terra::mask(bathy_raster, south_area) |>
  terra::trim() |>
  terra::disagg(fact = 10, method = "bilinear")

bathy_north <- terra::mask(bathy_raster, north_area) |>
  terra::trim() |>
  terra::disagg(fact = 10, method = "bilinear")

alaska_label <-
  data.frame(x = -160, y = 61, label = "Alaska") |>
  akgfmaps::transform_data_frame_crs(out.crs = "EPSG:3338")

p_main_map <- 
    ggplot() +
    geom_spatraster(data = bathy_raster, mapping = aes(fill = Height), alpha = 0.9) +
    geom_sf(data = map_layers$akland, color = NA, fill = "grey70") +
    geom_shadowtext(
      data = alaska_label,
                    mapping = aes(x = x, y = y, label = label), 
      color = "black", 
      bg.color = "white",
      size = 3.3
      ) +
    geom_sf(data = map_layers$survey.area, fill = NA, color = "black", linewidth = 0.5) +
    geom_shadowtext(
      data = area_label, 
      mapping = aes(x = x, y = y, label = SURVEY_ABBV), 
      color = "black", 
      bg.color = "white",
      size = 3.3) +
    geom_sf(data = station_loc, mapping = aes(color = factor(YEAR), shape = factor(YEAR))) +
    geom_sf(data = south_area, fill = NA, linetype = 2, color = "red", linewidth = 0.5) +
    geom_sf(data = north_area, fill = NA, linetype = 2, color = "red", linewidth = 0.5) +
    geom_sf(data = map_layers$graticule, alpha = 0.4, linewidth = 0.4, color = "grey30") +
    scale_x_continuous(
      limits = map_layers$plot.boundary$x + c(-1e5, 5e4),
                       breaks = map_layers$lon.breaks, 
      expand = c(0,0),
      oob = scales::oob_keep) +
    scale_y_continuous(
      limits = map_layers$plot.boundary$y + c(-1e5, 5e4),
                       breaks = map_layers$lat.breaks, 
      expand = c(0,0),
      oob = scales::oob_keep
      ) +
    scale_color_manual(name = "Year", values = c("2023" = "yellow", "2024" = "orange")) +
    scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17)) +
    scale_fill_distiller(
      name = "Depth (m)", 
      direction = 1, 
      breaks = c(0, 100, 200, 400, 600),
      labels = c("0", "100", "200", "400", "> 600"),
      limits = c(0, 600),
      oob = scales::squish) +
    theme_bw() +
    theme(axis.title = element_blank())

p_north <- 
    ggplot() +
    geom_spatraster(data = bathy_north, mapping = aes(fill = Height), alpha = 0.9) +
    geom_sf(data = map_layers$akland, color = NA, fill = "grey70") +
    geom_sf(data = map_layers$survey.area, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(
      data = dplyr::filter(station_loc, group == "North"), 
      mapping = aes(color = factor(YEAR), shape = factor(YEAR)),
      size = 1.5
    ) +
    geom_sf(data = map_layers$graticule, alpha = 0.4, linewidth = 0.3, color = "grey30") +
    scale_x_continuous(
      limits = c(bbox_north[['xmin']], bbox_north[['xmax']]) + c(-2.5e4 + 2.5e4),
      oob = scales::oob_keep,
      expand = c(0,0),
      breaks = c(-176, -177, -178)
    ) +
    scale_y_continuous(
      limits = c(bbox_north[['ymin']], bbox_north[['ymax']]) + c(-2.5e4 + 2.5e4),
      oob = scales::oob_keep,
      expand = c(0,0),
      breaks = c(59, 60)
    ) +
    scale_color_manual(name = "Year", values = c("2023" = "yellow", "2024" = "orange")) +
    scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17)) +
    scale_fill_distiller(
      name = "Depth (m)", 
      direction = 1, 
      breaks = c(0, 100, 200, 400, 600),
      labels = c("0", "100", "200", "400", "> 600"),
      limits = c(0, 600),
      oob = scales::squish) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      legend.position = "none")

p_south <- 
    ggplot() +
    geom_spatraster(data = bathy_south, mapping = aes(fill = Height), alpha = 0.9) +
    geom_sf(data = map_layers$akland, color = NA, fill = "grey70") +
    geom_sf(data = map_layers$survey.area, fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(data = dplyr::filter(station_loc, group == "South"), 
            mapping = aes(color = factor(YEAR), shape = factor(YEAR)),
            size = 1.5
    ) +
    geom_sf(data = map_layers$graticule, alpha = 0.4, linewidth = 0.7, color = "grey30") +
    scale_x_continuous(
      limits = c(bbox_south[['xmin']], bbox_south[['xmax']]) + c(-2.5e4 + 2.5e4),
      oob = scales::oob_keep,
      expand = c(0,0),
      breaks = c(-166, -167)
    ) +
    scale_y_continuous(
      limits = c(bbox_south[['ymin']], bbox_south[['ymax']]) + c(-2.5e4 + 2.5e4),
      oob = scales::oob_keep,
      expand = c(0,0),
      breaks = c(54.4, 54.8, 55.2)
    ) +
    scale_color_manual(name = "Year", values = c("2023" = "yellow", "2024" = "orange")) +
    scale_shape_manual(name = "Year", values = c("2023" = 16, "2024" = 17)) +
    scale_fill_distiller(
      name = "Depth (m)", 
      direction = 1, 
      breaks = c(0, 100, 200, 400, 600),
      labels = c("0", "100", "200", "400", "> 600"),
      limits = c(0, 600),
      oob = scales::squish) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      legend.position = "none")

p_map_grid <- 
  cowplot::plot_grid(
    p_main_map + 
      theme(
        axis.text = element_text(size = 9),
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9)
      ),
    cowplot::plot_grid(
      p_north + 
        theme(
          axis.text = element_text(size = 7.5),
          axis.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)
        ),
      p_south + 
        theme(
          axis.text = element_text(size = 7.5),
          axis.title = element_blank(),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)
        ),
      labels = c("B", "C"),
      nrow = 2
    ),
    labels = c("A", NULL),
    ncol = 2,
    rel_widths = c(0.6,0.35)
  )

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "map_shelf_slope.png"), 
    width = 169, 
    height = round(169/2.37), 
    units = "mm", 
    res = 300)
print(p_map_grid)
dev.off()
