library(akgfmaps)
library(ggspatial)
library(ggrepel)

slope_layers <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "WGS84")
shelf_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "WGS84")


project_areas <- sf::st_read(here::here("output", "2024_sample_zones_shelf_slope.shp"))
slope_allocation <- sf::st_read(here::here("analysis", "shelf_slope", "output", "2024_slope_allocation.shp"))

area_1 <- dplyr::filter(project_areas, label %in% c("Stratum 50", "Subarea 1"))
area_1_stations <- dplyr::filter(area_1, !is.na(STATIONID)) |>
  sf::st_centroid()
area_1_bbox <- area_1 |> sf::st_buffer(dist = 5000) |>
  sf::st_bbox()

allocation <- read.csv(file = here::here("analysis", "shelf_slope", "output", "slope_allocation.csv"))

allocation_points <- allocation |>
  sf::st_as_sf(coords = c("MEAN_LONGITUDE", "MEAN_LATITUDE"),
               crs = "WGS84")

dutch_harbor <- data.frame(x = -166.5417847, y = 53.8940888, label = "Dutch Harbor")


subarea1 <- ggplot() +
  geom_sf(data = slope_layers$akland) +
  geom_sf(data = shelf_layers$survey.grid, fill = NA) +
  geom_sf(data = slope_layers$survey.strata, fill = NA) +
  geom_sf(data = area_1,
          mapping = aes(fill = label),
          alpha = 0.5) +
  geom_sf_text(data = area_1_stations,
               mapping = aes(label = STATIONID)) +
  geom_sf(data = allocation_points,
          mapping = aes(shape = PRIMARY, color = PRIMARY)) +
  geom_text_repel(data = allocation,
            mapping = aes(x = MEAN_LONGITUDE, y = MEAN_LATITUDE, label = STATIONID)) +
  annotation_scale() +
  scale_fill_brewer(name = "Area") +
  scale_color_manual(name = "Priority", values = c("#000000", "#E69F00")) +
  ggtitle("Shelf/Slope Project Sample Areas\nSubarea 1 and Stratum 50") +
  scale_shape(name = "Priority") +
  scale_x_continuous(limits = c(area_1_bbox['xmin'], area_1_bbox['xmax'])) +
  scale_y_continuous(limits = c(area_1_bbox['ymin'], area_1_bbox['ymax'])) +
  theme_bw() +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9, 0.95))

ragg::agg_png(filename = here::here("plots"))
print(subarea1)

