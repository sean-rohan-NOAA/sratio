# Find new survey stations

library(sratio)
library(RODBC)
library(akgfmaps)
library(plotly)

con <- sratio::get_connected()

spp <- c(21740, 21720, 10210, 10261, 10110, 10112, 10130, 10285, 471, 68560, 68580, 69322)

map_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = "EPSG:3338")

stn_centroids <- map_layers$survey.grid |>
  sf::st_centroid()

stn_centroids <- sf::st_coordinates(stn_centroids) |>
  as.data.frame() |>
  dplyr::bind_cols(sf::st_drop_geometry(stn_centroids))

catch_2024 <- RODBC::sqlQuery(channel = con,
                              query = "select c.SPECIES_CODE,
                              c.VESSEL,
                              c.CRUISE,
                              c.HAUL,
                              c.REGION,
                              c.SPECIES_NAME,
                              c.TOTAL_WEIGHT,
                              c.TOTAL_NUMBERS,
                              h.STATION as STATIONID, 
                              h.HAUL_TYPE, 
                              h.PERFORMANCE, 
                              h.DISTANCE_FISHED_1 as DISTANCE_FISHED, 
                              h.GEAR_TEMPERATURE,
                              h.SURFACE_TEMPERATURE,
                              h.PERFORMANCE
                              from race_data.v_extract_edit_catch c,
                              race_data.v_extract_edit_haul h
                              where h.VESSEL = c.VESSEL
                              and h.CRUISE = c.CRUISE
                              and h.HAUL = c.HAUL
                              and h.HAUL_TYPE = 3
                              and h.CRUISE = 202401
                              and h.REGION = 'BS'
                              and h.PERFORMANCE >=0
                              ") |>
  dplyr::filter(SPECIES_CODE %in% spp)

hauls_2024 <- RODBC::sqlQuery(channel = con,
                              query = "select STATION as STATIONID 
                              from race_data.v_extract_edit_haul 
                              where HAUL_TYPE = 3
                              and CRUISE = 202401
                              and REGION = 'BS'
                              and PERFORMANCE >=0") |>
  dplyr::inner_join(stn_centroids)

sampled_stations <- sf::st_as_sf(hauls_2024, crs = "EPSG:3338", coords = c("X", "Y"))

spread_2024 <- read.csv(file = "C:\\Users\\sean.rohan\\Work\\afsc\\trawlmetrics\\output\\race_data_edit_hauls_table_EBS_2024.csv") |>
  dplyr::select(VESSEL, CRUISE, HAUL, EDIT_NET_SPREAD)

catch_2024 <- dplyr::inner_join(catch_2024, spread_2024) |> 
  dplyr::mutate(AREA_SWEPT_KM2 = DISTANCE_FISHED * EDIT_NET_SPREAD/1000) |>
  dplyr::mutate(CPUE_KG_KM2 = TOTAL_WEIGHT/AREA_SWEPT_KM2)

idw_results <- vector(mode = "list",
                      length = length(spp))

# Load sample zones
sample_zones <- readxl::read_xlsx(here::here("analysis", "15_30", "data", "expanded_sampling_zones.xlsx"))

sample_zones <- dplyr::inner_join(map_layers$survey.grid, sample_zones) |>
  dplyr::group_by(CLUSTER, TARGET) |>
  summarise(do_union = TRUE)

sample_zone_centroid <- sf::st_centroid(sample_zones)

pdf(file = here::here("analysis", "15_30", "plots", "map_idw_stations_2024.pdf"), 
    width = 8,
    height = 6,
    onefile = TRUE)

for(ii in 1:length(spp)) {
  
  sel_catch <- dplyr::filter(catch_2024, SPECIES_CODE == spp[ii]) |>
  dplyr::full_join(hauls_2024) |>
  dplyr::mutate(TOTAL_NUMBERS = if_else(is.na(TOTAL_NUMBERS), 0, TOTAL_NUMBERS),
                CPUE_KG_KM2 = if_else(is.na(CPUE_KG_KM2), 0, CPUE_KG_KM2),
                TOTAL_WEIGHT = if_else(is.na(TOTAL_WEIGHT), 0, TOTAL_WEIGHT)) |>
    dplyr::mutate(COMMON_NAME = sratio::species_code_label(spp[ii], type = "common_name"),
                  CPUE_KGHA = CPUE_KG_KM2) |>
    dplyr::rename(LONGITUDE = X,
                  LATITUDE = Y)
  
  idw_results[[ii]] <- akgfmaps::make_idw_map(x = sel_catch,
                         region = "sebs",
                         in.crs = "EPSG:3338",
                         key.title.units = "CPUE (kg/km2)",
                         extrapolation.grid.type = "sf")
  
  species_map <- idw_results[[ii]]$plot + 
    geom_sf(data = sampled_stations, fill = NA, color = "black", size = 0.4) +
      scale_x_continuous(limits = c(idw_results[[ii]]$map_layers$plot.boundary$x)) +
      scale_y_continuous(limits = c(idw_results[[ii]]$map_layers$plot.boundary$y))

  print(species_map)
}

dev.off()


# Load sample zones 

pdf(file = here::here("analysis", "15_30", "plots", "map_additional_zones_2024.pdf"), 
    width = 8,
    height = 6,
    onefile = TRUE)

for(ii in 1:length(spp)) {
  
  zone_map <- idw_results[[ii]]$plot +
    geom_sf(data = sampled_stations, fill = NA, color = "black", size = rel(1.2)) +
    geom_sf(data = sample_zones,
            mapping = aes(color = factor(CLUSTER)),
            fill = NA,
            alpha = 0.5,
            linewidth = 1.2) +
    geom_sf_text(data = sample_zone_centroid,
                 mapping = aes(label = CLUSTER,
                               color = factor(CLUSTER)),
                 face = "bold",
                 size = 5) +
    scale_color_tableau(name = "Cluster", guide = "none") +
    scale_x_continuous(limits = c(idw_results[[ii]]$map_layers$plot.boundary$x)) +
    scale_y_continuous(limits = c(idw_results[[ii]]$map_layers$plot.boundary$y))
  
  print(zone_map)
  
}

dev.off()

dir.create(here::here("analysis", "15_30", "output", "nav_files"), showWarnings = FALSE)

sf::st_write(obj = sample_zones, 
             dsn = here::here("analysis", "15_30", "output", "nav_files", "1530_bonus_zones.shp"))


# Sample zone polygons

sample_zones <- readxl::read_xlsx(here::here("analysis", "15_30", "data", "expanded_sampling_zones.xlsx"))