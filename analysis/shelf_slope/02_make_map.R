library(akgfmaps)
library(shadowtext)
library(cowplot)

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

study_area |>
  dplyr::rename(label = strat_label) |>
  sf::st_transform(crs = "WGS84") |>
  sf::st_write(here::here("output", "2024_sample_zones_shelf_slope.shp"), append = FALSE)

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


 # "#1AE4B6FF" "#FABA39FF" "#7A0403FF"

plot_strat <- function(x, slope_layers, sel_stratum, border_color = "#30123BFF") {
  
  x <- dplyr::filter(study_area, strat_label == sel_stratum)
  
  x_lab <- sf::st_centroid(x)
  
  focal_bbox <- sf::st_bbox(sf::st_buffer(x, dist = 5e4))
  
  broader_bbox <- sf::st_bbox(study_area)
  
  focal_plot <- ggplot() +
    geom_sf(data = slope_layers$bathymetry, color = "grey70", linewidth = 0.2) +
    geom_sf(data = slope_layers$akland,
            fill = "grey30", color = NA) +
    geom_sf(data = x, fill = NA, color = border_color) +
    geom_sf_text(data = x_lab, mapping = aes(label = STATIONID)) +
    scale_x_continuous(breaks = slope_layers$lon.breaks,
                       limits = c(focal_bbox['xmin'], focal_bbox['xmax'])) +
    scale_y_continuous(breaks = slope_layers$lat.breaks,
                       limits = c(focal_bbox['ymin'], focal_bbox['ymax'])) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.85),
          axis.title = element_blank(),
          legend.title = element_blank())
  
  strat_lab <- dplyr::group_by(x, strat_label) |>
    dplyr::summarise(do_union = TRUE) |>
    sf::st_centroid()
  
  strat_lab <- cbind(strat_lab, as.data.frame(sf::st_coordinates(strat_lab)))
  
  broad_plot <- ggplot() +
    geom_sf(data = slope_layers$bathymetry, color = "grey70", linewidth = 0.2) +
    geom_sf(data = slope_layers$akland,
            fill = "grey30", color = NA) +
    geom_sf(data = study_area, fill = NA, color = "grey50", alpha = 0.7) +
    geom_sf(data = x,
            fill = border_color, color = border_color) +
    geom_shadowtext(data = strat_lab, 
                       mapping = aes(x = X, y = Y, label = strat_label)) +
    scale_x_continuous(breaks = slope_layers$lon.breaks,
                       limits = c(broader_bbox['xmin'], broader_bbox['xmax'])) +
    scale_y_continuous(breaks = slope_layers$lat.breaks,
                       limits = c(broader_bbox['ymin'], broader_bbox['ymax'])) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.85),
          axis.title = element_blank(),
          legend.title = element_blank())
  
  return(list(focal_plot = focal_plot, broad_plot = broad_plot))

}



stratum_61 <- plot_strat(study_area, slope_layers = slope_layers, sel_stratum = "Stratum 61", border_color = "#30123BFF")

ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "focal_map_stratum61.png"), res = 180,
              width = 3, height = 3*1.29, units = "in")
print(stratum_61$focal_plot)
dev.off()


stratum_50 <- plot_strat(study_area, slope_layers = slope_layers, sel_stratum = "Stratum 50", border_color = "#30123BFF")

ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "focal_map_stratum50.png"), res = 180,
              width = 3, height = 3, units = "in")
print(stratum_50$focal_plot)
dev.off()
