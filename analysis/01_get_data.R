# Get data for selectivity ratio analysis

get_data <- function(species_codes, min_sample_size) {
  
  channel <- get_connected(schema = "AFSC")
  
  # Get haul data ----
  
  # 1995 hauls (manual selection). Not added into all_hauls.rds
  hauls_1995 <- RODBC::sqlQuery(channel = channel,
                                query = "select h.hauljoin, h.net_measured, h.wire_length, 
                            h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.region, h.duration, h.distance_fished, h.net_width, 
                            h.net_height, h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.surface_temperature, h.gear_temperature, h.haul_type 
                            from racebase.haul h  
                            where cruise = 199501
                            and vessel in (88, 89)
                            and haul_type = 7
                            and performance >= 0")|>
    dplyr::filter(
      (VESSEL == 88 & HAUL %in% c(192:196, 204, 205, 209, 211, 212, 215, 217, 223, 225, 228, 229, 234, 235, 237, 239, 245:246, 258:261))|
        (VESSEL == 89 & HAUL %in% c(207:211, 217, 218, 222, 224, 225, 228, 230:232, 234, 236, 239, 240, 246, 247, 251:252, 262, 263))) |>  #removed V89 H 253, 260; V88 H 247, 254 (bc out of station grid)
    #remove digits prior to station identification and format station numbers into double digits
    dplyr:: mutate(STATIONID = stringr::str_replace(gsub("^\\d+", "", STATIONID), "\\d+", function(x) sprintf('%02d', as.numeric(x)))) |>
    #hard code in new station ids for stations with multiple pairs to assign unique matchup id
    dplyr::mutate(STATIONID = 
                    case_when(
                      STATIONID == "A-06" & HAUL == 260 | HAUL == 261 ~ "A-06.a",
                      STATIONID == "A-06" & HAUL == 259 | HAUL == 258 ~ "A-06.b",
                      STATIONID == "B-08" & HAUL == 262 ~ "A-06.c",
                      STATIONID == "B-08" & HAUL == 263 ~ "A-06.c",
                      STATIONID == "B-08" & HAUL == 252 ~ "B-08.a",
                      STATIONID == "B-08" & HAUL == 246 & VESSEL == 88 ~ "B-08.a",
                      STATIONID == "B-08" & HAUL == 245 | HAUL == 251 ~ "B-08.b",
                      STATIONID == "F-04" & HAUL == 232 | HAUL == 231 ~ "F-04.a",
                      STATIONID == "F-04" & HAUL == 234 | HAUL == 223 ~ "F-04.b",
                      STATIONID == "F-04" & HAUL == 225 | HAUL == 236 ~ "F-04.c",
                      STATIONID == "F-05" & HAUL == 240 | HAUL == 229 ~ "F-05.a",
                      STATIONID == "F-05" & HAUL == 239 ~ "F-05.b",
                      STATIONID == "F-05" & HAUL == 228 & VESSEL == 88 ~ "F-05.b",
                      STATIONID == "F-05" & HAUL == 235 | HAUL == 234 ~ "F-05.c",
                      STATIONID == "F-08" & HAUL == 237 ~ "F-08.a", 
                      STATIONID == "F-08" & HAUL == 246 & VESSEL == 89 ~ "F-08.a",
                      STATIONID == "F-08" & HAUL == 247 ~ "F-08.b",
                      STATIONID == "G-08" & HAUL == 239 ~ "F-08.b",
                      STATIONID == "G-02" & HAUL == 212 | HAUL == 225 ~ "G-02.a",
                      STATIONID == "G-02" & HAUL == 215 | HAUL == 228 ~ "G-02.b",
                      STATIONID == "G-02" & HAUL == 217 | HAUL == 230 ~ "G-02.c",
                      STATIONID == "K-18" & HAUL == 218 ~ "K-18.a",
                      STATIONID == "J-18" & HAUL == 205 ~ "K-18.a",
                      STATIONID == "K-18" & HAUL == 217 ~ "K-18.b",
                      STATIONID == "J-18" & HAUL == 204 ~ "K-18.b",
                      STATIONID == "K-18" & HAUL == 224 ~ "K-18.c",
                      STATIONID == "J-18" & HAUL == 211 ~ "K-18.c",
                      STATIONID == "K-18" & HAUL == 222 ~ "K-18.d",
                      STATIONID == "J-18" & HAUL == 209 ~ "K-18.d",
                      STATIONID == "K-18" & HAUL == 196 | HAUL == 211 ~ "K-18.e",
                      STATIONID == "L-18" & HAUL == 195 | HAUL == 210 ~ "L-18.a",
                      STATIONID == "L-18" & HAUL == 194 | HAUL == 209 ~ "L-18.b",
                      STATIONID == "L-18" & HAUL == 192 | HAUL == 207 ~ "L-18.c",
                      STATIONID == "L-18" & HAUL == 193 | HAUL == 208 ~ "L-18.d",
                      TRUE ~ STATIONID
                    ))
  
  # 1998: Hauls in Bristol Bay had haul_type = 18, 15 minute hauls outside of Bristol Bay also had haul_type 18; need to include gear code because underbag experiments were also conducted during the same year using haul_type 18 
  special_tows_1998 <- RODBC::sqlQuery(channel = channel, query = "select hauljoin, net_measured, wire_length, start_time, performance, 
    vessel, cruise, haul, duration, distance_fished, 
    net_width, start_latitude, end_latitude, start_longitude, 
    end_longitude, stationid, gear_depth, bottom_depth, gear, 
    accessories, haul_type 
    from racebase.haul 
    where vessel in (88, 89) 
    and cruise = 199801
    and haul_type = 18
    and performance >= 0
    and gear = 44
    and accessories = 15")
  
  # Standard tows from 
  standard_tows_1998 <- RODBC::sqlQuery(channel = channel, query = "select hauljoin, net_measured, wire_length, start_time, performance, 
    vessel, cruise, haul, duration, distance_fished, 
    net_width, start_latitude, end_latitude, start_longitude, 
    end_longitude, stationid, gear_depth, bottom_depth, gear, 
    accessories, haul_type 
    from racebase.haul 
    where vessel in (88, 89) 
    and cruise = 199801
    and haul_type = 3
    and performance >= 0
    and gear = 44
    and accessories = 15") |>
    dplyr::inner_join(anti_join(otto_key_1998, 
                                special_tows_1998, 
                                by = c("VESSEL", "CRUISE", "HAUL")) |>
                        dplyr::select(VESSEL, CRUISE, HAUL, TOW_PAIR), 
                      by = c("VESSEL", "CRUISE", "HAUL"))
  
  hauls_1998 <- special_tows_1998 |>
    dplyr::inner_join(
      dplyr::select(otto_key_1998, VESSEL, CRUISE, HAUL, TOW_PAIR), 
      by = c("VESSEL", "CRUISE", "HAUL")
    ) |>
    dplyr::bind_rows(standard_tows_1998)
  
  hauls_1998$TOW_SPEED_KNOTS <- hauls_1998$DISTANCE_FISHED/hauls_1998$DURATION * 0.539957
  
  # Remove tows where tow speeds are outside of the acceptable range (2.8-3.2 knots)
  hauls_1998 <- dplyr::filter(hauls_1998, TOW_SPEED_KNOTS >= 2.8, TOW_SPEED_KNOTS <= 3.2)
  
  # Remove unpaired hauls
  hauls_1998 <- hauls_1998 |>
    dplyr::inner_join(
      as.data.frame(table(hauls_1998$VESSEL, hauls_1998$TOW_PAIR)) |>
        dplyr::rename(VESSEL = Var1, TOW_PAIR = Var2) |>
        dplyr::mutate(VESSEL = as.numeric(as.character(VESSEL)),
                      TOW_PAIR = as.numeric(as.character(TOW_PAIR))) |>
        dplyr::filter(Freq > 1) |>
        dplyr::select(-Freq),
      by = c("VESSEL", "TOW_PAIR")
    ) |>
    dplyr::mutate(STATIONID = paste0(VESSEL, "-", TOW_PAIR))
  
  # 2021: All hauls in 2021 had haul_type = 20
  hauls_2021 <- RODBC::sqlQuery(channel = channel,
                                query = "select h.hauljoin, h.net_measured, h.wire_length, h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.duration, h.distance_fished, h.net_width, 
                            h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.gear_temperature, h.haul_type 
                            from racebase.haul h 
                            where h.haul_type = 20
                            and h.performance >= 0
                            and h.cruise > 202100
                            and h.cruise < 202199
                            and h.region = 'BS'") |>
    dplyr::arrange(VESSEL, START_TIME)
  
  # 2021: Only stations with successful side-by-side comparisons
  hauls_2021 <- hauls_2021 |>
    dplyr::group_by(STATIONID) |>
    dplyr::summarise(n = n()) |>
    dplyr::filter(n > 1) |>
    dplyr::select(STATIONID) |>
    dplyr::inner_join(hauls_2021, 
                      by = c("STATIONID"))
  
  # 2022: 15 minute hauls were conducted alongside normal hauls in 2022 so comparison hauls have a mix of haul_types 3 and 20. 
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
                            and h.cruise > 202200
                            and h.cruise < 202299
                            and h.region = 'BS'") |>
    dplyr::arrange(VESSEL, START_TIME) |>
    dplyr::inner_join(unique(dplyr::select(hauls_short_2022, STATIONID, CRUISE)),
                      by = c("CRUISE", "STATIONID"))
  
  # 2023: 15 minute hauls were conducted opportunistically in 2023 at crab special project stations (L. Zacher). 
  # Haul types are a mix of 4 (15 minute) and 3
  hauls_short_2023 <- RODBC::sqlQuery(channel = channel,
                                       query = "select h.hauljoin, h.net_measured, h.wire_length, h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.region, h.duration, h.distance_fished, h.net_width, 
                            h.net_height, h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.surface_temperature, h.gear_temperature, h.haul_type 
                            from racebase.haul h 
                            where h.haul_type = 4
                            and performance >= 0
                            and h.cruise = 202301
                            and h.region = 'BS'") |>
    dplyr::arrange(VESSEL, START_TIME)
  
  hauls_normal_2023 <- RODBC::sqlQuery(channel = channel,
                                       query = "select h.hauljoin, h.net_measured, h.wire_length, h.start_time, 
                            h.performance, h.vessel, h.cruise, 
                            h.haul, h.region, h.duration, h.distance_fished, h.net_width, 
                            h.net_height, h.start_latitude, h.end_latitude, h.start_longitude, 
                            h.end_longitude, h.stationid, h.gear_depth, h.bottom_depth, h.gear, 
                            h.accessories, h.surface_temperature, h.gear_temperature, h.haul_type 
                            from racebase.haul h 
                            where h.haul_type = 3
                            and performance >= 0
                            and h.cruise = 202301
                            and h.region = 'BS'") |>
    dplyr::inner_join(unique(dplyr::select(hauls_short_2023, STATIONID, CRUISE)),
                      by = c("CRUISE", "STATIONID"))
    
  all_hauls <- dplyr::bind_rows(hauls_2021,
                                hauls_short_2022,
                                hauls_normal_2022,
                                hauls_short_2023,
                                hauls_normal_2023,
                                hauls_1998,
                                hauls_1995
                                )  |>
    dplyr::mutate(AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED) |>
    dplyr::mutate(TREATMENT = factor(
      plyr::round_any(DURATION, 0.25), 
      levels = c(0.25, 0.5), 
      labels = c(15, 30)
    )
    )
  
  no_pair <- dplyr::group_by(all_hauls, STATIONID, CRUISE, TREATMENT) |>
    dplyr::summarise(n = n(), .groups = "keep") |>
    dplyr::filter(n > 1) |>
    dplyr::select(-n)
  
  message("The following stations/cruises have no pair and will be removed:", no_pair)
  
  all_hauls <- dplyr::anti_join(all_hauls, no_pair, by = c("STATIONID", "CRUISE", "TREATMENT"))
  
  all_hauls <- dplyr::select(all_hauls, STATIONID, CRUISE) |>
    unique() |>
    dplyr::mutate(MATCHUP = dplyr::row_number()) |>
    dplyr::inner_join(all_hauls, by = c("STATIONID", "CRUISE")) |>
    dplyr::mutate(YEAR = floor(CRUISE/100))
  
  all_hauls |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n()) |>
    dplyr::select(YEAR = CRUISE, N_HAULS = n) |>
    dplyr::mutate(YEAR = floor(YEAR/100)) |>
    write.csv(here::here("plots", "n_hauls.csv"), row.names = FALSE)
  
  # Get catch data ----
  catch <- RODBC::sqlQuery(channel = channel,
                           query = paste0("select * from racebase.catch
                         where cruise > 199500
                         and region = 'BS'
                         and species_code in (", paste(species_codes, collapse = ","),  ")")) |>
    dplyr::inner_join(unique(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL)),
                      by = c("VESSEL", "CRUISE", "HAUL")) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, WEIGHT, NUMBER_FISH, HAULJOIN) |>
    dplyr::inner_join(dplyr::select(all_hauls, HAULJOIN, MATCHUP), by = "HAULJOIN") 
  
  # Get length data ----
  lengths <- RODBC::sqlQuery(channel = channel,
                             query = paste0("select * from racebase.length
                         where cruise in (", paste(use_cruises, collapse = ", "), 
                                            ") and region = 'BS'
                         and species_code in (", paste(species_codes, collapse = ","),  ")")) |>
    dplyr::inner_join(unique(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL)), 
                      by = c("VESSEL", "CRUISE", "HAUL")) |>
    dplyr::mutate(LENGTH = LENGTH/10) |>
    dplyr::select(VESSEL, CRUISE, HAUL, SPECIES_CODE, LENGTH, FREQUENCY, SEX, HAULJOIN) |>
    dplyr::inner_join(dplyr::select(all_hauls, HAULJOIN, MATCHUP), by = "HAULJOIN")
  
  # Check number of hauls
  dplyr::select(catch, VESSEL, CRUISE, HAUL) |>
    unique() |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n(), .groups = "keep")
  
  # Get crab carapace data ----
  # 1998: Get Goddard (1997) and Somerton et al. (2002) data
  crab_1995_1998 <- sratio::crab_size_1995_1998 |>
    dplyr::select(-STATIONID) |>
    dplyr::filter(!(SPECIES_CODE %in% c(69322, 69400))) |> # Remove BKC and horsehair crab
    dplyr::mutate(FREQUENCY = SAMPLING_FACTOR)
  
  
  # 2021 and 2022 Crab
  crab_2021_2022 <- RODBC::sqlQuery(channel = channel,
                                  paste0("select species_code, sex, shell_condition, length, width, weight, vessel, cruise, haul, sampling_factor frequency from crab.ebscrab_15_30_comparison_project
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
    dplyr::inner_join(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL, HAULJOIN) |>
                        unique(),
                      by = c("VESSEL", "CRUISE", "HAUL"))
  
  # 2023 Crab
  crab_2023 <- RODBC::sqlQuery(channel = channel,
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
    dplyr::inner_join(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL, HAULJOIN) |>
                        unique(),
                      by = c("VESSEL", "CRUISE", "HAUL"))
  
  crab <- dplyr::bind_rows(crab_2021_2022, crab_2023, crab_1995_1998) |>
    dplyr::inner_join(dplyr::select(all_hauls, HAULJOIN, MATCHUP), by = "HAULJOIN")
  
  crab$FREQUENCY[is.na(crab$FREQUENCY)] <- 1

  crab_fish <- dplyr::bind_rows(crab, lengths)
  
  # Identify matchups to use for selectivity analysis based on minimum sample size ----
  selectivity_flag <- dplyr::select(crab_fish, VESSEL, CRUISE, HAUL, MATCHUP, SPECIES_CODE, FREQUENCY) |>
    dplyr::group_by(VESSEL, CRUISE, HAUL, MATCHUP, SPECIES_CODE) |>
    dplyr::summarise(N_MEASURED = sum(FREQUENCY), .groups = "keep") |>
    dplyr::mutate(USE_FOR_SELECTIVITY = N_MEASURED > min_sample_size) |>
    dplyr::inner_join(dplyr::select(all_hauls, VESSEL, CRUISE, HAUL, MATCHUP, TREATMENT),
                      by = c("VESSEL", "CRUISE", "HAUL", "MATCHUP")) |>
    dplyr::ungroup() |>
    dplyr::group_by(MATCHUP, SPECIES_CODE) |>
    dplyr::summarise(USE_FOR_SELECTIVITY = sum(USE_FOR_SELECTIVITY), .groups = "keep") |>
    dplyr::mutate(USE_FOR_SELECTIVITY = USE_FOR_SELECTIVITY > 1)
  
  catch <- dplyr::inner_join(catch, 
                             selectivity_flag, 
                             by = c("SPECIES_CODE", "MATCHUP"))
  
  crab_fish <- dplyr::inner_join(crab_fish,
                                 selectivity_flag,
                                 by = c("SPECIES_CODE", "MATCHUP"))
  
  # test <- dplyr::select(catch, SPECIES_CODE, MATCHUP, USE_FOR_SELECTIVITY) |>
  #   unique() |>
  #   arrange(SPECIES_CODE, MATCHUP)
  # View(test)
  
  # Check number of hauls with size data
  dplyr::select(crab_fish, VESSEL, CRUISE, HAUL) |>
    unique() |>
    dplyr::group_by(CRUISE) |>
    dplyr::summarise(n = n(), .groups = "keep")
  
  data_1530 <- list(
    project = "15/30 Minute Tow Comparison",
    catch = catch,
    haul = all_hauls,
    size = crab_fish)

  save(data_1530, file = here::here("data", "data_1530.rda"))

  # saveRDS(object = all_hauls, file = here::here("data", "hauls_1530.rds"))
  # 
  # saveRDS(object = catch, file = here::here("data", "catch_1530.rds"))

  # saveRDS(object = crab_fish, file = here::here("data", "fish_crab_size_1530.rds"))
  
  crab_fish |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::filter(USE_FOR_SELECTIVITY) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY), .groups = "keep") |>
    write.csv(file = here::here("plots", "sample_sizes_1530.csv"), row.names = FALSE)
  
  crab_fish |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY), .groups = "keep") |>
    write.csv(file = here::here("plots", "sample_sizes_no_filter_1530.csv"), row.names = FALSE)

}

get_data(species_codes = c(21740, 21720, 10210, 10261, 10110, 10130, 10285, 471, 68560, 68580, 69322),
         min_sample_size = min_sample_size)
