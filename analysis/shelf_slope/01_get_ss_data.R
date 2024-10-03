library(sratio)

get_ss_data <- function(species_codes) {
  
  channel <- sratio::get_connected(schema = "AFSC")
  
  ss_haul_2023 <- RODBC::sqlQuery(channel = channel, 
                                  query = "select * from racebase.haul 
                                        where cruise = 202301 
                                        and haul_type = 20
                                        and vessel in (134, 162)
                                        and performance >= 0")  |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
                  TREATMENT = factor(GEAR, levels = c(172, 44)))
  
  ss_haul_akk_2024 <- RODBC::sqlQuery(channel = channel, 
                                  query = "select * from racebase.haul 
                                        where cruise = 202401 
                                        and haul_type = 20
                                        and vessel = 162
                                        and gear = 172
                                        and (haul < 18 or haul > 210)
                                        and accessories = 115
                                        and performance >= 0")  |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
                  TREATMENT = factor(GEAR, levels = c(172, 44)))
  
  ss_haul_nwe_2024 <- RODBC::sqlQuery(channel = channel, 
                                      query = "select * from racebase.haul 
                                        where cruise = 202401 
                                        and vessel = 134
                                        and performance >= 0")  |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
                  TREATMENT = factor(GEAR, levels = c(172, 44))) |>
    dplyr::filter(STATIONID %in% ss_haul_akk_2024$STATIONID)
  
  ss_haul_akk_2024 <- dplyr::filter(ss_haul_akk_2024, 
                                    STATIONID %in% ss_haul_nwe_2024$STATIONID)
  
  ss_hauls <- dplyr::bind_rows(ss_haul_2023, ss_haul_akk_2024, ss_haul_nwe_2024)
  
  
  ss_matchups <- dplyr::group_by(ss_hauls, STATIONID, CRUISE) |>
    dplyr::summarise(n = n()) |>
    dplyr::filter(n > 1) |>
    dplyr::select(-n) |>
    unique() |>
    dplyr::ungroup() |>
    dplyr::mutate(MATCHUP = dplyr::row_number())
  
  ss_hauls <- dplyr::inner_join(ss_hauls, ss_matchups) |>
    dplyr::mutate(YEAR = floor(CRUISE/100))
  
  ss_catch <- RODBC::sqlQuery(channel = channel, 
                                   query = paste0("select * from racebase.catch 
                                        where cruise in (202301, 202401) 
                                        and vessel in (134, 162)
                                        and species_code in (", paste(species_codes, collapse = ","), ")",
                                        " and hauljoin in (", paste(ss_hauls$HAULJOIN, collapse = ","), ")")) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, WEIGHT, NUMBER_FISH, HAULJOIN) |>
    dplyr::inner_join(dplyr::select(ss_hauls, MATCHUP, HAULJOIN) |>
                        unique())
  
  ss_length <- RODBC::sqlQuery(channel = channel, 
                                    query = paste0("select * from racebase.length 
                                        where cruise in (202301, 202401)
                                        and vessel in (134, 162)
                                        and species_code in (", paste(species_codes, collapse = ","), ")",
                                                   " and hauljoin in (", paste(ss_hauls$HAULJOIN, collapse = ","), ")")) |>
    dplyr::mutate(LENGTH = LENGTH/10) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, LENGTH, FREQUENCY, SEX, HAULJOIN) |>
    dplyr::inner_join(dplyr::select(ss_hauls, MATCHUP, HAULJOIN) |>
                        unique())
  
  ss_crab_2023 <- RODBC::sqlQuery(channel = channel,
                               paste0("select species_code, sex, shell_condition, length, width, weight, vessel, cruise, haul, sampling_factor frequency from crab.ebscrab_15_30_slope_shelf_comparison_2023tows
                  where species_code in (", paste(species_codes, collapse = ","), ")")) |>
    dplyr::mutate(VESSEL = as.numeric(VESSEL),
                  CRUISE = as.numeric(CRUISE),
                  HAUL = as.numeric(HAUL),
                  SPECIES_CODE = as.numeric(SPECIES_CODE),
                  SEX = as.numeric(SEX),
                  SHELL_CONDITION = as.numeric(SHELL_CONDITION),
                  LENGTH = as.numeric(LENGTH),
                  WIDTH = as.numeric(WIDTH),
                  FREQUENCY = as.numeric(FREQUENCY)) |>
    dplyr::inner_join(dplyr::select(ss_hauls, VESSEL, CRUISE, HAUL, MATCHUP, HAULJOIN) |>
                        unique())
  
  ss_crab_2024 <- read.csv(file = here::here("analysis", "shelf_slope", "data", "ebscrab_15_30_slope_shelf_comparison_2024tows.csv")) |>
    dplyr::select(-WEIGHT) |>
    dplyr::inner_join(dplyr::select(ss_hauls, VESSEL, CRUISE, HAUL, MATCHUP, HAULJOIN) |>
                        unique())
  
  ss_crab_2024 <- ss_crab_2024[, which(names(ss_crab_2024) %in% names(ss_crab_2023))]
  
  ss_crab <- dplyr::bind_rows(ss_crab_2023, ss_crab_2024)
  
  ss_crab$FREQUENCY[is.na(ss_crab$FREQUENCY)] <- 1
  
  ss_crab_fish <- dplyr::bind_rows(ss_crab, ss_length) |>
    dplyr::filter(HAULJOIN %in% ss_hauls$HAULJOIN,
                  SPECIES_CODE %in% species_codes)
  
  data_ss <- list(
    project = "Shelf/Slope Tow Comparison",
    catch = ss_catch,
    haul = ss_hauls,
    size = ss_crab_fish)
  
  save(data_ss, file = here::here("data", "data_ss.rda"))
  
  saveRDS(ss_hauls, file = here::here("analysis", "shelf_slope", "data", "ss_haul.rds"))
  saveRDS(ss_catch, file = here::here("analysis", "shelf_slope", "data", "ss_catch.rds"))
  saveRDS(ss_length, file = here::here("analysis", "shelf_slope", "data", "ss_length.rds"))
  saveRDS(object = ss_crab, file = here::here("analysis", "shelf_slope", "data", "ss_crab.rds"))
  saveRDS(object = ss_crab_fish, file = here::here("analysis", "shelf_slope", "data", "ss_fish_crab_size.rds"))
  
  ss_crab_fish |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY)) |>
    write.csv(file = here::here("analysis", "shelf_slope", "plots", "ss_sample_sizes.csv"), row.names = FALSE)
  
  ss_hauls |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n()) |>
    dplyr::select(YEAR = CRUISE, N_HAULS = n) |>
    dplyr::mutate(YEAR = floor(YEAR/100)) |>
    write.csv(here::here("analysis", "shelf_slope", "plots", "n_hauls.csv"), row.names = FALSE)
  
}

get_ss_data(species_codes = c(21740, 21720, 10130, 10115, 10110, 10112, 471, 68580, 658560, 69322))
