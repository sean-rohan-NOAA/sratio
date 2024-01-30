library(sratio)
library(ggthemes)

channel <- sratio::get_connected(schema = "AFSC")

n_stratum_61 <- 16
n_stratum_11 <- 12

bad_substrate <- function(performance) {
  
  results <- sapply(abs(performance), 
                    function(x) (x > 0 & x < 3) | (x > 3.99 & x < 4.4))
  
  return(results)
  
}

slope_lut <- readxl::read_xlsx(here::here("analysis", "shelf_slope", "docs", "all_slope_stations new.xlsx")) |>
  dplyr::filter(CRUISE < 200800) |>
  dplyr::select(VESSEL, CRUISE, HAUL, STRATUM = old_STRATUM, STATIONID = Old_STATIONID, NEW_STRATUM = new_strata, NEW_STATIONID = New_StationID)

slope_hauls <- RODBC::sqlQuery(channel = channel, 
                               query = "select rbh.cruise, rbh.vessel, rbh.haul, rbh.gear_depth, rbh.bottom_depth, rbh.start_time, 
                               rbh.start_latitude, rbh.start_longitude, rbh.end_latitude, rbh.end_longitude, 
                               rbh.performance, rbh.haul_type, rbh.stratum, rbh.stationid 
                               from race_data.cruises c, race_data.surveys s, racebase.haul rbh 
                               where rbh.cruise = c.cruise 
                               and s.survey_id = c.survey_id 
                               and c.vessel_id = rbh.vessel 
                               and s.survey_definition_id = 78 
                               and rbh.stratum in (11, 61)")

slope_before_2008 <- dplyr::filter(slope_hauls, CRUISE < 200800) |>
  dplyr::inner_join(slope_lut) |>
  dplyr::mutate(STRATUM = as.numeric(NEW_STRATUM),
                STATIONID = NEW_STATIONID) |>
  dplyr::select(-NEW_STRATUM, -NEW_STATIONID)

slope_after_2008 <- dplyr::filter(slope_hauls, CRUISE > 200800)

slope_hauls <- dplyr::bind_rows(slope_before_2008, slope_after_2008) |>
  dplyr::mutate(BAD_SUBSTRATE = bad_substrate(PERFORMANCE))


slope_summary <- dplyr::group_by(slope_hauls, STATIONID, STRATUM) |>
  dplyr::summarise(MEAN_BOTTOM_DEPTH = mean(BOTTOM_DEPTH), 
                   MEAN_LATITUDE = mean(START_LATITUDE),
                   MEAN_LONGITUDE = mean(START_LONGITUDE),
                   n_total = n(),
                   n_bad = sum(BAD_SUBSTRATE)) |>
  dplyr::mutate(p_bad = n_bad/n_total) |>
  dplyr::arrange(STRATUM, -p_bad) |>
  dplyr::filter(n_total >= 2, p_bad == 0)

bssa1_layers <- akgfmaps::get_base_layers(select.region = "bssa1", set.crs = "WGS84")
bssa6_layers <- akgfmaps::get_base_layers(select.region = "bssa6", set.crs = "WGS84")
slope_layers <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "WGS84")
shelf_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "WGS84")

slope_stations_sf <- sf::st_as_sf(slope_summary, 
                                  coords = c("MEAN_LONGITUDE", "MEAN_LATITUDE"),
                                  crs = "WGS84")

set.seed(19740)

project_stations <- dplyr::bind_rows(
  data.frame(STATIONID = sample(slope_summary$STATIONID[slope_summary$STRATUM == 61], replace = FALSE), 
             STRATUM = 61,
             SAMPLE_PRIORITY = 1:length(unique(slope_summary$STATIONID[slope_summary$STRATUM == 61]))),
  data.frame(STATIONID = sample(slope_summary$STATIONID[slope_summary$STRATUM == 11], replace = FALSE), 
             STRATUM = 11,
             SAMPLE_PRIORITY = 1:length(unique(slope_summary$STATIONID[slope_summary$STRATUM == 11])))
) |>
  dplyr::inner_join(slope_summary)

project_stations_sf$PRIMARY <- NA

project_stations$PRIMARY[project_stations$STRATUM == 61] <- ifelse(project_stations$SAMPLE_PRIORITY[project_stations$STRATUM == 61] <= n_stratum_61, "Primary", "Alternate")
project_stations$PRIMARY[project_stations$STRATUM == 11] <- ifelse(project_stations$SAMPLE_PRIORITY[project_stations$STRATUM == 11] <= n_stratum_11, "Primary", "Alternate")

project_stations_sf <- dplyr::inner_join(slope_stations_sf, project_stations)

p1 <- ggplot() +
  geom_sf(data = shelf_layers$bathymetry) +
  geom_sf(data = slope_layers$survey.strata, alpha = 0.5) +
  geom_sf(data = slope_layers$akland) +
  geom_sf(data = project_stations_sf,
          mapping = aes(color = PRIMARY,
                        shape = factor(STRATUM))) +
  scale_x_continuous(limits = slope_layers$plot.boundary$x) +
  scale_y_continuous(limits = slope_layers$plot.boundary$y) +  
  scale_shape(name = "Stratum") +
  scale_color_colorblind(name = "Priority") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.25))

plot_stratum_11 <- ggplot() +
  geom_sf(data = shelf_layers$bathymetry) +
  geom_sf(data = slope_layers$akland) +
  geom_sf(data = bssa1_layers$survey.strata, alpha = 0.5) +
  geom_sf(data = dplyr::filter(project_stations_sf, STRATUM == 11),
          mapping = aes(fill = PRIMARY),
          shape = 21) +
  scale_fill_colorblind(name = "Priority") +
  scale_x_continuous(limits = bssa1_layers$plot.boundary$x) +
  scale_y_continuous(limits = bssa1_layers$plot.boundary$y) +  
  theme_bw() +
  theme(legend.position = c(0.85, 0.85))

plot_stratum_61 <- ggplot() +
  geom_sf(data = shelf_layers$bathymetry) +
  geom_sf(data = slope_layers$akland) +
  geom_sf(data = bssa6_layers$survey.strata, alpha = 0.5) +
  geom_sf(data = dplyr::filter(project_stations_sf, STRATUM == 61),
          mapping = aes(fill = PRIMARY),
          shape = 21) +
  scale_fill_colorblind(name = "Priority") +
  scale_x_continuous(limits = bssa6_layers$plot.boundary$x) +
  scale_y_continuous(limits = bssa6_layers$plot.boundary$y) +  
  theme_bw() +
  theme(legend.position = c(0.15, 0.12))

png(filename = here::here("analysis", "shelf_slope", "plots", "sampling_2024_stations.png"), 
              width = 5, height = 6, units = "in", res = 300)
print(p1)
dev.off()

png(filename = here::here("analysis", "shelf_slope", "plots", "sampling_2024_stations_subarea1.png"), 
    width = 5.6, height = 4, units = "in", res = 300)
print(plot_stratum_11)
dev.off()

png(filename = here::here("analysis", "shelf_slope", "plots", "sampling_2024_stations_subarea6.png"), 
    width = 4.2, height = 5, units = "in", res = 300)
print(plot_stratum_61)
dev.off()

write.csv(x = project_stations, file = here::here("analysis", "shelf_slope", "output", "slope_allocation.csv"), row.names = FALSE)
