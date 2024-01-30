library(akgfmaps)

slope_layers <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "EPSG:3338")

shelf_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")

stn_pattern <- c(paste0(LETTERS[12:18], "-31"), 
                 paste0(LETTERS[12:18], "-32"), 
                 "AZ0504", "A-04", "A-03", "A-02", "B-03", "B-02", "B-01", "C-18", "C-01")

study_area <- shelf_layers$survey.grid[grepl(paste(stn_pattern, collapse="|"), shelf_layers$survey.grid$STATIONID), ] |>
  sf::st_join(shelf_layers$survey.strata) |>
  sf::st_make_valid() |>
  dplyr::mutate(strat_label = paste0("Stratum ", Stratum),
                STRATUm = Stratum) |>
  dplyr::bind_rows(dplyr::filter(slope_layers$survey.strata, STRATUM %in% c(61, 11)) |>
                     dplyr::group_by(STRATUM) |>
                     dplyr::summarise(do_union = TRUE) |>
                     dplyr::mutate(strat_label = paste0("Subarea ", floor(STRATUM/10)))
                   ) |>
  dplyr::select(strat_label, STRATUM, STATIONID)

study_bbox <- sf::st_bbox(study_area)

ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "ebs_slope_strata.png"), width = 120, height = 120 * 1.21,
              units = "mm", res = 300)
print(
ggplot() +
  geom_sf(data = slope_layers$akland,
          fill = "grey30") +
  geom_sf(data = slope_layers$bathymetry, color = "grey40") +
  geom_sf(data = slope_layers$survey.strata |>
            dplyr::mutate(subarea = factor(floor(STRATUM/10))),
          mapping = aes(fill = subarea),
          color = NA,
          alpha = 0.8) +
  coord_sf(xlim = slope_layers$plot.boundary$x,
           ylim = slope_layers$plot.boundary$y) +
  scale_x_continuous(breaks = slope_layers$lon.breaks) +
  scale_y_continuous(breaks = slope_layers$lat.breaks) +
  scale_fill_viridis_d(name = "Subarea", option = "H") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
)
dev.off()


ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "shelf_slope_target_areas.png"), width = 120, height = 120 * 1.21,
              units = "mm", res = 300)
print(
ggplot() +
  geom_sf(data = slope_layers$akland,
          fill = "grey30", color = NA) +
  geom_sf(data = study_area,
          mapping = aes(fill = strat_label), color = "black", alpha = 0.7) +
  # geom_sf_text(data = sf::st_centroid(study_area),
  #              mapping = aes(label = STATIONID)) +
  geom_sf(data = slope_layers$bathymetry, color = "grey40", linewidth = 0.2) +

  scale_x_continuous(breaks = slope_layers$lon.breaks,
                     limits = c(study_bbox['xmin'], study_bbox['xmax'])) +
  scale_y_continuous(breaks = slope_layers$lat.breaks,
                     limits = c(study_bbox['ymin'], study_bbox['ymax'])) +
  scale_fill_viridis_d(name = "Subarea", option = "H") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.85),
        axis.title = element_blank(),
        legend.title = element_blank())
)
dev.off()

