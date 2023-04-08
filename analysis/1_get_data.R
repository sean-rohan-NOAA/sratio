# Get data for selectivity ratio analysis

get_data <- function(species_codes) {
  
  channel <- get_connected(schema = "AFSC")
  
  # Get haul data ----
  
  # All hauls in 2021 had haul_type = 20
  hauls_2021 <- RODBC::sqlQuery(channel = channel,
                                query = "select h.hauljoin, h.net_measured, h.wire_length, h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.region, h.duration, h.distance_fished, h.net_width, 
                            h.net_height, h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.surface_temperature, h.gear_temperature, h.haul_type 
                            from racebase.haul h 
                            where h.haul_type = 20
                            and h.performance >= 0
                            and h.cruise > 202100
                            and h.cruise < 202199
                            and h.region = 'BS'") |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED) |>
    dplyr::arrange(VESSEL, START_TIME)
  
  # Only stations with successful side-by-side comparisons
  hauls_2021 <- hauls_2021 |>
    dplyr::group_by(STATIONID) |>
    dplyr::summarise(n = n()) |>
    dplyr::filter(n > 1) |>
    dplyr::select(STATIONID) |>
    dplyr::inner_join(hauls_2021)
  
  # 15 minute hauls were conducted alongside normal hauls in 2022 so comparison hauls have a mix of haul_types 3 and 20. 
  hauls_short_2022 <- RODBC::sqlQuery(channel = channel,
                                      query = "select h.hauljoin, h.net_measured, h.wire_length, h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.region, h.duration, h.distance_fished, h.net_width, 
                            h.net_height, h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.surface_temperature, h.gear_temperature, h.haul_type 
                            from racebase.haul h 
                            where h.haul_type = 20
                            and h.cruise > 202200
                            and h.cruise < 202299
                            and h.region = 'BS'") |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED) |>
    dplyr::arrange(VESSEL, START_TIME)
  
  hauls_normal_2022 <- RODBC::sqlQuery(channel = channel,
                                       query = "select h.hauljoin, h.net_measured, h.wire_length, h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.region, h.duration, h.distance_fished, h.net_width, 
                            h.net_height, h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.surface_temperature, h.gear_temperature, h.haul_type 
                            from racebase.haul h 
                            where h.haul_type = 3
                            and performance >= 0
                            and h.cruise > 202100
                            and h.region = 'BS'") |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED) |>
    dplyr::arrange(VESSEL, START_TIME) |>
    dplyr::inner_join(unique(dplyr::select(hauls_short_2022, STATIONID, CRUISE)))
  
  all_hauls <- dplyr::bind_rows(hauls_2021,
                                hauls_short_2022,
                                hauls_normal_2022) |>
    dplyr::mutate(TREATMENT = factor(
      plyr::round_any(DURATION, 0.25), 
      levels = c(0.25, 0.5), 
      labels = c(15, 30)
    )
    )
  
  
  no_pair <- dplyr::group_by(all_hauls, STATIONID, CRUISE, TREATMENT) |>
    dplyr::summarise(n = n()) |>
    dplyr::filter(n > 1) |>
    dplyr::select(-n)
  
  message("The following stations/cruises have no pair and will be removed:", no_pair)
  
  
  all_hauls <- dplyr::anti_join(all_hauls, no_pair)
  
  all_hauls <- dplyr::select(all_hauls, STATIONID, CRUISE) |>
    unique() |>
    dplyr::mutate(MATCHUP = dplyr::row_number()) |>
    dplyr::inner_join(all_hauls)
  
  saveRDS(object = all_hauls, file = here::here("data", "hauls_1530.rds"))
  
  all_hauls |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n()) |>
    dplyr::select(YEAR = CRUISE, N_HAULS = n) |>
    dplyr::mutate(YEAR = floor(YEAR/100)) |>
    write.csv(here::here("plots", "n_hauls.csv"), row.names = FALSE)
  
  
  # Get catch data ----
  catch <- RODBC::sqlQuery(channel = channel,
                           query = paste0("select * from racebase.catch
                         where cruise > 202100
                         and region = 'BS'
                         and species_code in (", paste(species_codes, collapse = ","),  ")")) |>
    dplyr::inner_join(unique(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL))) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, WEIGHT, NUMBER_FISH, HAULJOIN)
  
  saveRDS(object = catch, file = here::here("data", "catch_1530.rds"))
  
  # Check number of hauls
  dplyr::select(catch, VESSEL, CRUISE, HAUL) |>
    unique() |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n())
  
  # Get length data ----
  lengths <- RODBC::sqlQuery(channel = channel,
                             query = paste0("select * from racebase.length
                         where cruise > 202100
                         and region = 'BS'
                         and species_code in (", paste(species_codes, collapse = ","),  ")")) |>
    dplyr::inner_join(unique(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL))) |>
    dplyr::mutate(LENGTH = LENGTH/10) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, LENGTH, FREQUENCY, SEX, HAULJOIN)
  
  saveRDS(object = lengths, file = here::here("data", "fish_lengths_1530.rds"))
  
  # Check number of hauls with length data
  dplyr::select(lengths, VESSEL, CRUISE, HAUL) |>
    unique() |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n())
  
  
  lengths |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY)) |>
    write.csv(file = here::here("plots", "sample_sizes_1530.csv"), row.names = FALSE)
  
  
  # Get crab carapace data ----
  
}

get_data(species_codes = species_codes)
