library(sratio)
library(ggthemes)
library(ggrepel)

channel <- sratio::get_connected(schema = "AFSC")

n_stratum_61 <- c("61-17", "61-04", "61-03", "61-20", "61-22", "61-06", "61-05", "61-07", "61-21", "61-12", "61-18", "61-13", "61-14", "61-15")
n_stratum_11 <- c("11-21", "11-20", "11-05", "11-06", "11-23", "11-07", "11-08", "11-36", "11-25", 
                  "11-33", "11-27", "11-37", "11-28", "11-13")

skate_hapc <- c("61-02", "11-18")

bad_substrate <- function(performance) {
  
  results <- sapply(abs(performance), 
                    function(x) (x > 0 & x < 3) | (x > 3.99 & x < 4.4))
  
  return(results)
  
}

slope_lut <- readxl::read_xlsx(here::here("analysis", "shelf_slope", "docs", "all_slope_stations_new.xlsx")) |>
  dplyr::filter(CRUISE < 200800) |>
  dplyr::select(VESSEL, CRUISE, HAUL, STRATUM = old_STRATUM, STATIONID = Old_STATIONID, NEW_STRATUM = new_strata, NEW_STATIONID = New_StationID)

slope_hauls <- RODBC::sqlQuery(channel = channel, 
                               query = "select rbh.cruise, rbh.vessel, rbh.haul, rbh.gear_depth, rbh.bottom_depth, rbh.start_time, 
                               rbh.start_latitude, rbh.start_longitude, rbh.end_latitude, rbh.end_longitude, 
                               rbh.performance, rbh.haul_type, rbh.stratum, rbh.stationid, hpn.note 
                               from race_data.cruises c, race_data.surveys s, racebase.haul rbh, race_data.haul_performance_codes hpc, race_data.haul_performance_notes hpn 
                               where rbh.cruise = c.cruise 
                               and s.survey_id = c.survey_id 
                               and c.vessel_id = rbh.vessel 
                               and s.survey_definition_id = 78
                               and hpn.haul_performance_note_id = hpc.haul_performance_note_id
                               and rbh.performance = hpc.haul_performance_code
                               and rbh.stratum in (11, 61)") |>
  dplyr::filter(!(STATIONID %in% skate_hapc))

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

project_stations$PRIMARY <- ifelse(project_stations$STATIONID %in% c(n_stratum_61, n_stratum_11) , "Primary", "Alternate")

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


project_station_labels <- sf::st_coordinates(project_stations_sf) |>
  dplyr::bind_cols(project_stations_sf) |>
  as.data.frame()


plot_stratum_11_labels <- ggplot() +
  geom_sf(data = shelf_layers$bathymetry) +
  geom_sf(data = slope_layers$akland) +
  geom_sf(data = dplyr::mutate(bssa1_layers$survey.strata, ID = 1) |>
            dplyr::group_by(ID) |>
            dplyr::summarise(do_union = TRUE), alpha = 0.5, color = "grey50", fill = "grey50") +
  geom_sf(data = dplyr::filter(project_stations_sf, STRATUM == 11),
          mapping = aes(fill = PRIMARY, shape = PRIMARY)) +
  geom_text_repel(data = dplyr::filter(project_station_labels, STRATUM == 11),
          mapping = aes(x = X, y = Y, label = STATIONID), size = rel(2)) +
  scale_fill_colorblind(name = "Priority") +
  scale_x_continuous(limits = bssa1_layers$plot.boundary$x) +
  scale_y_continuous(limits = bssa1_layers$plot.boundary$y) +  
  scale_shape_manual(name = "Priority", values = c(21,22)) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.82),
        axis.title = element_blank(),
        axis.text = element_text(size = 7.5))

plot_stratum_61_labels <- ggplot() +
  geom_sf(data = shelf_layers$bathymetry) +
  geom_sf(data = slope_layers$akland) +
  geom_sf(data = dplyr::mutate(bssa6_layers$survey.strata, ID = 1) |>
            dplyr::group_by(ID) |>
            dplyr::summarise(do_union = TRUE), alpha = 0.5, color = "grey50", fill = "grey50") +
  geom_sf(data = dplyr::filter(project_stations_sf, STRATUM == 61),
          mapping = aes(fill = PRIMARY, shape = PRIMARY)) +
  geom_text_repel(data = dplyr::filter(project_station_labels, STRATUM == 61),
                  mapping = aes(x = X, y = Y, label = STATIONID), size = rel(2)) +
  scale_fill_colorblind(name = "Priority") +
  scale_x_continuous(limits = bssa6_layers$plot.boundary$x) +
  scale_y_continuous(limits = bssa6_layers$plot.boundary$y) +  
  scale_shape_manual(name = "Priority", values = c(21,22)) +
  theme_bw() +
  theme(legend.position = c(0.80, 0.82),
        axis.title = element_blank())


# Make shapefiles for navigation software
project_stations_sf |>
  dplyr::select(STATIONID, STRATUM, DEPTH = MEAN_BOTTOM_DEPTH, PRIORITY = SAMPLE_PRIORITY, PRIMARY) |>
  sf::st_write(here::here("analysis", "shelf_slope", "output", "2024_slope_allocation.shp"), 
               append = FALSE)

slope_hauls |>
  sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"), crs = "WGS84") |>
  dplyr::select(VESSEL, CRUISE, HAUL, TIME = START_TIME, PERFORMANCE, PERFDESC = NOTE, STATIONID) |>
  sf::st_write(here::here("analysis", "shelf_slope", "output", "slope_tow_starts.shp"), 
               append = FALSE)

dplyr::bind_rows(slope_hauls |>
                   sf::st_as_sf(coords = c("START_LONGITUDE", "START_LATITUDE"), crs = "WGS84"),
                 slope_hauls |>
                   sf::st_as_sf(coords = c("END_LONGITUDE", "END_LATITUDE"), crs = "WGS84")
                 ) |>
  dplyr::select(VESSEL, CRUISE, HAUL, TIME = START_TIME, PERFORMANCE, PERFDESC = NOTE, STATIONID) |> 
  dplyr::group_by(VESSEL, CRUISE, HAUL, PERFORMANCE, PERFDESC, STATIONID, TIME) |> 
  dplyr::summarize(do_union = FALSE) |> 
  sf::st_cast("LINESTRING") |>
  sf::st_write(here::here("analysis", "shelf_slope", "output", "slope_towpaths.shp"), 
               append = FALSE)

# Make plots
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


png(filename = here::here("analysis", "shelf_slope", "plots", "sampling_2024_stations_subarea1_lab.png"), 
    width = 4.2, height = 3, units = "in", res = 300)
print(plot_stratum_11_labels)
dev.off()

png(filename = here::here("analysis", "shelf_slope", "plots", "sampling_2024_stations_subarea6_lab.png"), 
    width = 3.2, height = 4, units = "in", res = 300)
print(plot_stratum_61_labels)
dev.off()

write.csv(x = project_stations, 
          file = here::here("analysis", "shelf_slope", "output", "slope_allocation.csv"), 
          row.names = FALSE)
