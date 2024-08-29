library(akgfmaps)
library(shadowtext)
library(cowplot)
library(sratio)

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

ss_tow_locations <- sratio::data_ss$haul |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"),
               crs = "WGS84") |>
  sf::st_transform(crs = "EPSG:3338") |>
  dplyr::mutate(label = "Shelf/Slope\nComparison Tow")

plot_strat_samples <- function(x, slope_layers, sel_stratum, tow_locations, legend_position) {
  
  x <- dplyr::filter(study_area, strat_label == sel_stratum)
  
  x_lab <- sf::st_centroid(x)
  
  focal_bbox <- sf::st_bbox(sf::st_buffer(x, dist = 5e4))
  
  broader_bbox <- sf::st_bbox(study_area)
  
  focal_plot <- ggplot() +
    geom_sf(data = slope_layers$bathymetry, 
            color = "grey50", 
            linewidth = 0.2) +
    geom_sf(data = slope_layers$survey.strata, 
            fill = NA, 
            color = "grey20", 
            linewidth = 0.3) +
    geom_sf(data = shelf_layers$survey.strata, 
            fill = NA, 
            color = "grey20", 
            linewidth = 0.3) +
    geom_sf(data = slope_layers$akland,
            color = NA) +
    geom_sf(data = tow_locations,
            mapping = aes(color = label),
            size = 2) + 
    scale_x_continuous(breaks = slope_layers$lon.breaks,
                       limits = c(focal_bbox['xmin'], focal_bbox['xmax'])) +
    scale_y_continuous(breaks = slope_layers$lat.breaks,
                       limits = c(focal_bbox['ymin'], focal_bbox['ymax'])) +
    scale_color_manual(values = "#309BD3") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white", colour = "grey20"), 
          strip.background = element_rect(fill = "grey85",colour = "grey20"), 
          legend.spacing.y = unit(-0.35, "cm"),
          legend.text = element_text(size = 7),
          legend.background=element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = legend_position,
          legend.box.just = "left",
          legend.box = "vertical", 
          axis.text = element_text(size = 8),
          legend.title = element_blank(),
          axis.title = element_blank())
  
  return(focal_plot)
  
}


subarea_6 <- plot_strat_samples(study_area, 
                                slope_layers = slope_layers, 
                                sel_stratum = "Stratum 61",
                                tow_locations = ss_tow_locations,
                                legend_position = c(0.25, 0.08))

ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "subarea6_samples.png"), res = 180,
              width = 3, height = 3*1.29, units = "in")
print(subarea_6)
dev.off()


subarea_1 <- plot_strat_samples(study_area, 
                        slope_layers = slope_layers, 
                        sel_stratum = "Stratum 50",
                        tow_locations = ss_tow_locations,
                        legend_position = c(0.25, 0.08))

ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "subarea1_samples.png"), res = 180,
              width = 3, height = 3, units = "in")
print(subarea_1)
dev.off()


ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "subarea1_samples_multi.png"), res = 180,
              width = 6, height = 3, units = "in")
print(
  cowplot::plot_grid(
    subarea_1 + 
      ggtitle("Subarea 1/Stratum 50") + 
      theme(title = element_text(hjust = 0,
                                 size = 12)),
    subarea_6 + 
      ggtitle("Subarea 6/Stratum 61") + 
      theme(title = element_text(hjust = 0,
                                 size = 12),
            legend.position = "none"),
    align = "v"
  )
)
dev.off()
