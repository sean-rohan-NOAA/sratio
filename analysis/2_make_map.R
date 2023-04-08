# Make plot of stations
library(akgfmaps)

haul_loc_df <- readRDS(file = here::here("data", "hauls_1530.rds")) |>
  dplyr::mutate(LATITUDE = (START_LATITUDE + END_LATITUDE) /2,
                LONGITUDE = (START_LONGITUDE + END_LONGITUDE) /2,
                YEAR = floor(CRUISE/100)) |>
  dplyr::select(STATIONID, YEAR, VESSEL, BOTTOM_DEPTH, LATITUDE, LONGITUDE) |>
  dplyr::group_by(STATIONID, YEAR) |>
  dplyr::summarise(LATITUDE = mean(LATITUDE),
                   LONGITUDE = mean(LONGITUDE),
                   BOTTOM_DEPTH = mean(BOTTOM_DEPTH)) |>
  sf::st_as_sf(crs = "EPSG:4326", coords = c("LONGITUDE", "LATITUDE")) |>
  sf::st_transform(crs = "EPSG:3338")


map_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

ragg::agg_png(filename = here::here("plots", "sample_map.png"), width = 4, height = 3, units = "in", res = 300)
print(
ggplot() +
  geom_sf(data = map_layers$akland, 
          linewidth = 0.2, 
          fill = "grey50", 
          color = "black") +
  geom_sf(data = map_layers$survey.area, 
          fill = NA, 
          color = "black") +
  geom_sf_text(data = sf::st_centroid(map_layers$survey.area),
               mapping = aes(label = factor(SURVEY, labels = c("EBS", "NBS")))) +
  geom_sf(data = map_layers$bathymetry, 
          color = "black",
          linewidth = 0.3) +
  geom_sf(data = haul_loc_df,
          mapping = aes(shape = factor(YEAR))) +
  geom_sf(data = map_layers$graticule, linewidth = 0.2, alpha = 0.7, color = "grey60") +
  coord_sf(xlim = map_layers$plot.boundary$x,
           ylim = map_layers$plot.boundary$y) +
  scale_x_continuous(breaks = map_layers$lon.breaks) +
  scale_y_continuous(breaks = map_layers$lat.breaks) +
  scale_shape(name = "Year", solid = FALSE) +
  theme_bw() +
  theme(axis.title = element_blank())
)
dev.off()
