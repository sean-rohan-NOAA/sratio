library(sratio)

get_ss_data <- function(species_codes) {
  
  channel <- sratio::get_connected(schema = "AFSC")
  
  ss_haul_2023 <- RODBC::sqlQuery(channel = channel, 
                                  query = "select * from racebase.haul 
                                        where cruise = 202301 
                                        and haul_type = 20
                                        and vessel in (134, 162)
                                        and performance >= 0")  |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED)
  
  ss_matchups_2023 <- dplyr::group_by(ss_haul_2023, STATIONID, CRUISE) |>
    dplyr::summarise(n = n()) |>
    dplyr::filter(n > 1) |>
    dplyr::select(-n) |>
    unique() |>
    dplyr::ungroup() |>
    dplyr::mutate(MATCHUP = dplyr::row_number())
  
  ss_haul_2023 <- dplyr::inner_join(ss_haul_2023, ss_matchups_2023)
  
  
  ss_matchups_2023 <- dplyr::select(ss_haul_2023, STATIONID, CRUISE) |>
    unique() |>
    dplyr::mutate(MATCHUP = dplyr::row_number()) |>
    dplyr::inner_join(
      dplyr::select(ss_haul_2023, VESSEL, CRUISE, HAUL)
    )
  
  ss_catch_2023 <- RODBC::sqlQuery(channel = channel, 
                                   query = paste0("select * from racebase.catch 
                                        where cruise = 202301 
                                        and vessel in (134, 162)
                                        and species_code in (", paste(species_codes, collapse = ","),  ")")) |>
    dplyr::inner_join(
      dplyr::select(ss_haul_2023, VESSEL, CRUISE, HAUL)
    ) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, WEIGHT, NUMBER_FISH, HAULJOIN)
  
  
  ss_length_2023 <- RODBC::sqlQuery(channel = channel, 
                                    query = paste0("select * from racebase.length 
                                        where cruise = 202301 
                                        and vessel in (134, 162)
                                        and species_code in (", paste(species_codes, collapse = ","),  ")")) |>
    dplyr::inner_join(
      dplyr::select(ss_haul_2023, VESSEL, CRUISE, HAUL)
    ) |>
    dplyr::mutate(LENGTH = LENGTH/10) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, LENGTH, FREQUENCY, SEX, HAULJOIN)
  
  saveRDS(ss_haul_2023, file = here::here("analysis", "shelf_slope", "data", "ss_haul.rds"))
  saveRDS(ss_catch_2023, file = here::here("analysis", "shelf_slope", "data", "ss_catch.rds"))
  saveRDS(ss_length_2023, file = here::here("analysis", "shelf_slope", "data", "ss_length.rds"))
  
  ss_length_2023 |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY)) |>
    write.csv(file = here::here("analysis", "shelf_slope", "plots", "sample_sizes_1530.csv"), row.names = FALSE)
  
  ss_haul_2023 |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n()) |>
    dplyr::select(YEAR = CRUISE, N_HAULS = n) |>
    dplyr::mutate(YEAR = floor(YEAR/100)) |>
    write.csv(here::here("analysis", "shelf_slope", "plots", "n_hauls.csv"), row.names = FALSE)
  
}

get_ss_data(species_codes = species_codes)
