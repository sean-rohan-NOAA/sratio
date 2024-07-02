# Make plot of stations
library(akgfmaps)
library(sratio)

# Load built-in data sets
haul_df <- sratio::data_1530$haul


haul_loc_df <- haul_df  |>
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

ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "map_samples_by_stratum.png"), 
              width = 6, height = 4, units = "in", res = 300)
print(
ggplot() +
  geom_sf(data = map_layers$akland, 
          linewidth = 0.2, 
          fill = "grey50", 
          color = "black") +
  geom_sf(data = map_layers$survey.strata, 
          fill = NA, 
          color = "black") +
  geom_sf_text(data = sf::st_centroid(map_layers$survey.strata),
               mapping = aes(label = Stratum)) +
  geom_sf(data = haul_loc_df,
          mapping = aes(color = factor(YEAR),
                        shape = factor(YEAR))) +
  geom_sf(data = map_layers$graticule, linewidth = 0.2, alpha = 0.7, color = "grey60") +
  coord_sf(xlim = map_layers$plot.boundary$x,
           ylim = map_layers$plot.boundary$y) +
  scale_x_continuous(breaks = map_layers$lon.breaks) +
  scale_y_continuous(breaks = map_layers$lat.breaks) +
  ggthemes::scale_color_colorblind(name = "Year") +
  scale_shape(name = "Year") +
  theme_bw() +
  theme(axis.title = element_blank())
)
dev.off()


# Multi-panel sample map
map_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")

haul_locs <- sratio::data_1530$haul |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"), crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338")


crab_layers <- dplyr::bind_rows(sf::st_read(here::here("assets", "Pribilof_RKC_strata.shp")) |>
                                  dplyr::mutate(area = "Pribilof Islands"),
                                sf::st_read(here::here("assets", "Pribilof_RKC_strata.shp")) |>
                                  dplyr::mutate(area = "St. Matthew  Island"),
                                sf::st_read(here::here("assets", "BBRKC_strata.shp")) |>
                                  dplyr::mutate(area = "Bristol Bay"),
                                sf::st_read(here::here("assets", "EBS_CO_CB_strata.shp")) |>
                                  dplyr::mutate(area = "EBS")
)


ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "sample_map_1995_2023.png"), 
              width = 169, height = 120, res = 300, units = "mm")
print(
  ggplot() +
    geom_sf(data = map_layers$akland) +
    geom_sf(data = map_layers$survey.strata, fill = NA) +
    geom_sf(data = haul_locs, size = rel(0.8), color = "#49C1ADFF") +
    facet_wrap(~YEAR) +
    scale_x_continuous(limits = map_layers$plot.boundary$x,
                       breaks = map_layers$lon.breaks) +
    scale_y_continuous(limits = map_layers$plot.boundary$y,
                       breaks = map_layers$lat.breaks) +
    theme_bw()
)
dev.off()


ragg::agg_png(filename = here::here("analysis", "15_30", "plots", "sample_map_crab_strata_1995_2023.png"), 
              width = 169, height = 120, res = 300, units = "mm")
print(
  ggplot() +
    geom_sf(data = map_layers$akland) +
    geom_sf(data = crab_layers, fill = NA) +
    geom_sf(data = haul_locs, size = rel(0.8), color = "#49C1ADFF") +
    facet_wrap(~YEAR) +
    scale_x_continuous(limits = map_layers$plot.boundary$x,
                       breaks = map_layers$lon.breaks) +
    scale_y_continuous(limits = map_layers$plot.boundary$y,
                       breaks = map_layers$lat.breaks) +
    theme_bw()
)
dev.off()
