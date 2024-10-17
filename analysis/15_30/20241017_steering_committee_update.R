library(sratio)
library(akgfmaps)
library(shadowtext)

# 15/30 comparison

sratio::data_1530$haul |>
  dplyr::select(YEAR, MATCHUP) |>
  unique() |>
  dplyr::group_by(YEAR) |>
  dplyr::summarise(n = n()) |>
  write.csv(file = here::here("analysis", "15_30", "plots",  "n_matchups.csv"), row.names = FALSE)


map_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")

hauls_1530 <- sratio::data_1530$haul |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"), crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338") |>
  dplyr::inner_join(data.frame(VESSEL = c(88, 89, 94, 134, 162),
                               VESSEL_NAME = c("Arcturus", "Aldebaran", "Vesteraalen", "Northwest Explorer", "Alaska Knight")))

sample_map <- ggplot() +
  geom_sf(data = map_layers$akland) +
  geom_sf(data = map_layers$survey.strata, linewidth = 0.1, alpha = 0.5, fill = NA) +
  geom_sf(data = hauls_1530,
          color = "deepskyblue2") +
  geom_sf(data = map_layers$graticule, 
          linewidth = 0.1,
          alpha = 0.3) +
  facet_wrap(~YEAR) +
  scale_x_continuous(limits = map_layers$plot.boundary$x,
                     breaks = map_layers$lon.breaks) +
  scale_y_continuous(limits = map_layers$plot.boundary$y,
                     breaks = map_layers$lat.breaks) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank())


ragg::agg_png(file = here::here("analysis", "15_30", "plots", "sample_map_by_year.png"), 
              width = 6, 
              height = 4, 
              units = "in",
              res = 300)
print(sample_map)
dev.off()



# Size measurements

sample_size_df <- read.csv(file = here::here("analysis", "15_30", "plots", "sample_sizes_wide_1530.csv")) |> 
  tidyr::pivot_longer(cols = 2:7) |>
  dplyr::rename(YEAR = name, N = value) |>
  dplyr::mutate(YEAR = as.numeric(gsub(pattern = "X", replacement = "", x = YEAR)))


# Shelf/slope comparison ----



# Setup analysis species
gear_names <- data.frame(GEAR = c(44, 172),
                         GEAR_NAME = c("83-112", "PNE"))

# Get base layers
bssa_layers <- akgfmaps::get_base_layers(select.region = "ebs.slope", 
                                         set.crs = "EPSG:3338")

bssa1_layers <- akgfmaps::get_base_layers(select.region = "bssa1", 
                                          set.crs = "EPSG:3338")

bssa6_layers <- akgfmaps::get_base_layers(select.region = "bssa6", 
                                          set.crs = "EPSG:3338")

shelf_layers <- akgfmaps::get_base_layers(select.region = "sebs",
                                          set.crs = "EPSG:3338")

target_strata <- bssa_layers$survey.strata |>
  dplyr::filter(STRATUM %in% c(11, 61)) |>
  dplyr::select(-B5_) |>
  dplyr::group_by(STRATUM) |>
  dplyr::summarise(do_union = TRUE) |>
  dplyr::mutate(STRATUM = paste("Slope Stratum ", STRATUM)) |>
  dplyr::bind_rows(
    shelf_layers$survey.strata |>
      dplyr::filter(Stratum == 61) |> 
      dplyr::select(STRATUM = Stratum) |>
      dplyr::mutate(STRATUM = paste("Shelf Stratum", STRATUM))
  )

survey_areas <- bssa_layers$survey.area |>
  dplyr::mutate(SURVEY = "Slope") |> 
  dplyr::bind_rows(dplyr::mutate(shelf_layers$survey.area, SURVEY = "Shelf"))

survey_areas$SURVEY <- factor(survey_areas$SURVEY, levels = c("Shelf", "Slope"))

haul_df <- sratio::data_ss$haul |>
  dplyr::inner_join(gear_names) |>
  dplyr::filter(CRUISE == 202401, GEAR == 44) |>
  dplyr::mutate(TYPE = 'Shelf/slope comparison tows') |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"),
               crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338")
  
  # Setup boundaries
  bssa1_x <- c(bssa1_layers$plot.boundary$x[1], bssa1_layers$plot.boundary$x[1] + 182000)
  bssa1_y <- c(bssa1_layers$plot.boundary$y[1], bssa1_layers$plot.boundary$y[1] + 183000)
  bssa6_x <- c(bssa6_layers$plot.boundary$x[1], bssa6_layers$plot.boundary$x[1] + 182000)
  bssa6_y <- c(bssa6_layers$plot.boundary$y[1], bssa6_layers$plot.boundary$y[1] + 183000) - 5000
  
  
  bathy <- dplyr::filter(bssa_layers$bathymetry, METERS %in% c(200, 400, 600, 800, 1000))
  
  bssa1_labels <- data.frame(x = c(-165.5, -165.9, -167.1, -167.6, -167.8), 
                             y = c(54.5, 54.55, 54.5, 54.5, 54.65),
                             label = c("200 m", "400 m", "600 m", "800 m", "1000 m")) |>
    akgfmaps::transform_data_frame_crs(out.crs = "EPSG:3338")
  
  bssa6_labels <- data.frame(x = c(-177.5, -178.4), 
                             y = c(59.25, 59.3),
                             label = c("200 m", "1000 m")) |>
    akgfmaps::transform_data_frame_crs(out.crs = "EPSG:3338")
  
  ss_map_subarea1 <- ggplot() +
    geom_sf(data = bssa_layers$akland) +
    geom_sf(data = bathy, linewidth = 0.1) +
    geom_sf(data = target_strata,
            mapping = aes(fill = STRATUM),
            alpha = 0.5,
            color = NA) +
    geom_sf(data = haul_df,
               mapping = aes(color = TYPE)) +
    geom_shadowtext(data = bssa1_labels,
                    mapping = aes(x = x, y = y, label = label),
                    size = 2.5,
                    bg.color = "white",
                    color = "black") +
    geom_sf(data = bssa6_layers$graticule, 
            alpha = 0.1, 
            linewidth = 0.3) +
    scale_color_manual(values = "black") +
    scale_fill_viridis_d(option = "E") +
    scale_x_continuous(limits = bssa1_x,
                       breaks = bssa1_layers$lon.breaks) +
    scale_y_continuous(limits = bssa1_y,
                       breaks = bssa1_layers$lat.breaks) +
    theme_bw() +
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          legend.direction = "horizontal")
  
  ss_map_subarea6 <- ggplot() +
    geom_sf(data = bssa_layers$akland) +
    geom_sf(data = bathy, linewidth = 0.1) +
    geom_sf(data = target_strata,
            mapping = aes(fill = STRATUM),
            alpha = 0.5,
            color = NA) +
    geom_sf(data = haul_df,
            mapping = aes(color = TYPE)) +
    geom_shadowtext(data = bssa6_labels,
                    mapping = aes(x = x, y = y, label = label),
                    size = 2.5,
                    bg.color = "white",
                    color = "black") +
    geom_sf(data = bssa6_layers$graticule, 
            alpha = 0.1, 
            linewidth = 0.3) +
    scale_color_manual(values = "black") +
    scale_fill_viridis_d(option = "E") +
    scale_x_continuous(limits = bssa6_x,
                       breaks = bssa6_layers$lon.breaks) +
    scale_y_continuous(limits = bssa6_y,
                       breaks = bssa6_layers$lat.breaks) +
    theme_bw() +
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank())
  
  plot_legend <- cowplot::get_legend(ss_map_subarea1 + theme(plot.margin = unit(c(0, 5, 2, 5), units = "mm")))
  
  
  ss_map_all_tows_2024 <- cowplot::plot_grid(
    cowplot::plot_grid(
      ss_map_subarea1 + theme(legend.position = "none", plot.margin = unit(c(2, 5, 9, 5), units = "mm")),
      ss_map_subarea6 + theme(legend.position = "none", plot.margin = unit(c(2, 5, 9, 5), units = "mm")),
      nrow = 1
      ),
      cowplot::plot_grid(NULL, 
                         plot_legend, 
                         NULL, ncol = 3, nrow = 1, rel_widths = c(0.1, 0.9, 0.1)),
      rel_heights = c(0.9, 0.1),
      nrow = 2)
  
  ragg::agg_png(here::here("analysis", "shelf_slope", "plots", "Shelf_slope_tow_locations_2024.png"),  
                height = 3, width = 6, units = "in", res = 300)
  print(ss_map_all_tows_2024)
  dev.off()
