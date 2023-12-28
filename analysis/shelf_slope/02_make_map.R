library(akgfmaps)

slope_layers <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "EPSG:3338")

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

