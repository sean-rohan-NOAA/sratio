# Get data for selectivity ratio analysis

get_data <- function(species_codes) {
  
  library(sratio)
  library(gapindex)
  
  channel <- sratio::get_connected(schema = "AFSC")
  
  # Get haul data ----
  
  # 1995 hauls (manual selection). Not added into all_hauls.rds
  hauls_1995 <- RODBC::sqlQuery(
    channel = channel, 
    query = "SELECT
              HAULJOIN,
              NET_MEASURED,
              WIRE_LENGTH,
              START_TIME,
              PERFORMANCE,
              VESSEL,
              CRUISE,
              HAUL,
              REGION,
              DURATION,
              DISTANCE_FISHED,
              NET_WIDTH,
              NET_HEIGHT,
              START_LATITUDE,
              START_LONGITUDE,
              END_LATITUDE,
              END_LONGITUDE,
              STATIONID,
              GEAR_DEPTH,
              BOTTOM_DEPTH,
              GEAR,
              ACCESSORIES,
              SURFACE_TEMPERATURE,
              GEAR_TEMPERATURE,
              HAUL_TYPE
            FROM 
              RACEBASE.HAUL H
            WHERE 
              CRUISE = 199501
              AND VESSEL IN (88, 89)
              AND HAUL_TYPE = 7
              AND PERFORMANCE >= 0")|>
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
  special_tows_1998 <- 
    RODBC::sqlQuery(
      channel = channel, 
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          REGION, 
          ACCESSORIES, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          VESSEL IN (88, 89) 
          AND CRUISE = 199801 
          AND HAUL_TYPE = 18 
          AND PERFORMANCE >= 0 
          AND GEAR = 44 
          AND ACCESSORIES = 15")
  
  # Standard tows from 
  standard_tows_1998 <- 
    RODBC::sqlQuery(
      channel = channel, 
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          REGION, 
          ACCESSORIES, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          VESSEL IN (88, 89) 
          AND CRUISE = 199801 
          AND HAUL_TYPE = 3 
          AND PERFORMANCE >= 0 
          AND GEAR = 44 
          AND ACCESSORIES = 15"
      ) |>
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
  hauls_2021 <- 
    RODBC::sqlQuery(
      channel = channel, 
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          SURFACE_TEMPERATURE, 
          REGION, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL H 
        WHERE 
          HAUL_TYPE = 20 
          AND PERFORMANCE >= 0 
          AND CRUISE > 202100 
          AND CRUISE < 202199 
          AND REGION = 'BS'"
      ) |>
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
  hauls_short_2022 <- 
    RODBC::sqlQuery(
      channel = channel,
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 20 
          AND CRUISE > 202200 
          AND CRUISE < 202299 
          AND REGION = 'BS'"
      ) |>
    dplyr::arrange(VESSEL, START_TIME)
  
  hauls_normal_2022 <- 
    RODBC::sqlQuery(
      channel = channel, 
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 3 
          AND PERFORMANCE >= 0 
          AND CRUISE > 202200 
          AND CRUISE < 202299 
          AND REGION = 'BS'"
      ) |>
    dplyr::arrange(VESSEL, START_TIME) |>
    dplyr::inner_join(unique(dplyr::select(hauls_short_2022, STATIONID, CRUISE)),
                      by = c("CRUISE", "STATIONID"))
  
  # 2023: 15 minute hauls were conducted opportunistically in 2023 at crab special project stations (L. Zacher). 
  # Haul types are a mix of 4 (15 minute) and 3
  hauls_short_2023 <- 
    RODBC::sqlQuery(
      channel = channel, 
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 4 
          AND PERFORMANCE >= 0 
          AND CRUISE = 202301 
          AND REGION = 'BS'"
      ) |>
    dplyr::arrange(VESSEL, START_TIME)
  
  hauls_normal_2023 <- 
    RODBC::sqlQuery(
      channel = channel,
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 3 
          AND PERFORMANCE >= 0 
          AND CRUISE = 202301 
          AND REGION = 'BS'"
      ) |>
    dplyr::inner_join(unique(dplyr::select(hauls_short_2023, STATIONID, CRUISE)),
                      by = c("CRUISE", "STATIONID"))
  
  # 2024: 15 minute hauls were conducted alongside index station sampling during legs 1-3 
  hauls_short_2024 <- 
    RODBC::sqlQuery(
      channel = channel,
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 20 
          AND GEAR = 44 
          AND HAUL < 210 
          AND PERFORMANCE >= 0 
          AND CRUISE = 202401 
          AND REGION = 'BS'"
      ) |>
    dplyr::arrange(VESSEL, START_TIME) |>
    dplyr::filter(STATIONID %in% akgfmaps::get_survey_stations(select.region = "sebs", include.corners = TRUE))
  
  hauls_normal_2024 <- 
    RODBC::sqlQuery(
      channel = channel,
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 3 
          AND GEAR = 44 
          AND PERFORMANCE >= 0 
          AND CRUISE = 202401 
          AND REGION = 'BS'"
      ) |>
    dplyr::inner_join(unique(dplyr::select(hauls_short_2024, STATIONID, CRUISE)),
                      by = c("CRUISE", "STATIONID"))
  
  bonus_hauls_2024 <- 
    RODBC::sqlQuery(
      channel = channel,
      query = 
        "SELECT 
          HAULJOIN, 
          NET_MEASURED, 
          WIRE_LENGTH, 
          START_TIME, 
          PERFORMANCE, 
          VESSEL, 
          CRUISE, 
          HAUL, 
          REGION, 
          DURATION, 
          DISTANCE_FISHED, 
          NET_WIDTH, 
          NET_HEIGHT, 
          START_LATITUDE, 
          END_LATITUDE, 
          START_LONGITUDE, 
          END_LONGITUDE, 
          STATIONID, 
          GEAR_DEPTH, 
          BOTTOM_DEPTH, 
          GEAR, 
          ACCESSORIES, 
          SURFACE_TEMPERATURE, 
          GEAR_TEMPERATURE, 
          HAUL_TYPE 
        FROM 
          RACEBASE.HAUL 
        WHERE 
          HAUL_TYPE = 20 
          AND HAUL > 210 
          AND GEAR = 44 
          AND PERFORMANCE >= 0 
          AND CRUISE = 202401 
          AND REGION = 'BS'"
      ) |>
    dplyr::filter(STATIONID %in% akgfmaps::get_survey_stations(select.region = "sebs", include.corners = TRUE))
  
  
  bonus_hauls_2024 <- 
    bonus_hauls_2024 |>
    dplyr::group_by(STATIONID) |>
    dplyr::summarise(N = n()) |>
    dplyr::filter(N > 1) |>
    dplyr::inner_join(bonus_hauls_2024, by = "STATIONID")
  
  all_hauls <- 
    dplyr::bind_rows(
      hauls_2021,
      hauls_short_2022,
      hauls_normal_2022,
      hauls_short_2023,
      hauls_normal_2023,
      hauls_normal_2024,
      hauls_short_2024,
      bonus_hauls_2024,
      hauls_1998,
      hauls_1995
    )  |>
    dplyr::mutate(
      AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
      TOW_SPEED_KNOTS = DISTANCE_FISHED/DURATION * 0.539957
    ) |>
    dplyr::mutate(
      TREATMENT = 
        factor(
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
    write.csv(here::here("analysis", "15_30", "plots", "n_hauls.csv"), row.names = FALSE)
  
  # Get net numbers ----
  
  net_number_1995_1998 <- RODBC::sqlQuery(
    channel = channel,
    query = paste0(
    "SELECT 
      HPDEN.NET_NUMBER, 
      HPDEN.VESSEL, 
      HPDEN.CRUISE, 
      HPDEN.HAUL
    FROM 
      RACE_EDIT.RB2_HPDEN HPDEN, 
      RACEBASE.HAUL RBH
    WHERE RBH.HAUL = HPDEN.HAUL
      AND RBH.VESSEL = HPDEN.VESSEL
      AND RBH.CRUISE = HPDEN.CRUISE
      AND RBH.HAULJOIN IN (", paste(all_hauls$HAULJOIN, collapse = ","), ")"
                                          )) |>
    dplyr::mutate(NET_NUMBER = as.numeric(NET_NUMBER),
                  VESSEL = as.numeric(VESSEL),
                  CRUISE = as.numeric(CRUISE),
                  HAUL = as.numeric(HAUL))
  
  net_number_2021_2024 <- RODBC::sqlQuery(
    channel = channel,
    query = paste0(
      "SELECT 
        RDH.NET_NUMBER,
        RDH.HAUL,
        RDC.VESSEL_ID VESSEL,
        RDC.CRUISE
      FROM 
        RACE_DATA.HAULS RDH,
        RACE_DATA.CRUISES RDC,
        RACEBASE.HAUL RBH
      WHERE
        RDC.CRUISE_ID = RDH.CRUISE_ID
        AND RDC.CRUISE = RBH.CRUISE
        AND RDH.HAUL = RBH.HAUL
        AND RDC.VESSEL_ID = RBH.VESSEL
        AND RBH.CRUISE >= 202100
        AND RBH.HAULJOIN IN (", paste(all_hauls$HAULJOIN, collapse = ","), ")"
    )
  )
  
  net_numbers <- dplyr::bind_rows(net_number_1995_1998, net_number_2021_2024)
  
  all_hauls <- dplyr::left_join(
    all_hauls,
    net_numbers,
    by = c("VESSEL", "CRUISE", "HAUL")
  )
  
  # Get catch data ----
  catch <- RODBC::sqlQuery(
    channel = channel,
    query = paste0(
      "SELECT HAULJOIN, 
        VESSEL, 
        CRUISE, 
        HAUL, 
        SPECIES_CODE, 
        WEIGHT, 
        NUMBER_FISH 
      FROM 
        RACEBASE.CATCH 
      WHERE 
        CRUISE > 199500
        AND REGION = 'BS'
        AND SPECIES_CODE IN (", paste(species_codes, collapse = ","),  ") 
        AND HAULJOIN IN (", paste(all_hauls$HAULJOIN, collapse = ","), ")")) |>
    dplyr::inner_join(all_hauls[c("HAULJOIN", "MATCHUP")], by = "HAULJOIN") 
  
  # Get length data ----
  fish_lengths <- 
    RODBC::sqlQuery(
    channel = channel,
    query = paste0(
      "SELECT
          VESSEL, 
          CRUISE,
          HAUL,
          SPECIES_CODE,
          LENGTH,
          FREQUENCY,
          SEX,
          HAULJOIN
        FROM 
          RACEBASE.LENGTH
        WHERE 
          REGION = 'BS'
          AND SPECIES_CODE IN (", paste(species_codes, collapse = ","),  ") 
          AND HAULJOIN IN (", paste(all_hauls$HAULJOIN, collapse = ","), ")"
      )
    ) |>
    dplyr::mutate(LENGTH = LENGTH/10) |>
    dplyr::inner_join(all_hauls[c("HAULJOIN", "MATCHUP")], by = "HAULJOIN")
  
  # Calculate sampling factors for fish
  fish_lengths <- 
    fish_lengths |>
    dplyr::group_by(HAULJOIN, SPECIES_CODE) |>
    dplyr::summarise(N_LENGTHS = sum(FREQUENCY)) |>
    dplyr::inner_join(catch[c("HAULJOIN", "SPECIES_CODE", "NUMBER_FISH")], 
                      by = c("SPECIES_CODE", "HAULJOIN")) |>
    dplyr::ungroup() |>
    dplyr::mutate(SAMPLING_FACTOR = NUMBER_FISH/N_LENGTHS) |>
    dplyr::select(HAULJOIN, SPECIES_CODE, SAMPLING_FACTOR) |>
    dplyr::inner_join(fish_lengths) |>
    dplyr::rename(SIZE = LENGTH)
  
  fish_lengths$SAMPLING_FACTOR[fish_lengths$SAMPLING_FACTOR < 1] <- 1
  
  # Get crab carapace data ----
  # Provided by Shannon Hennessey on April 14, 2025
  # Includes data from some hauls in 1995 that were not included because they were not included in a tow pair
  crab <- 
    read.csv(file = here::here("analysis", "15_30", "data", "specimen_1530_ss.csv")) |>
    dplyr::filter(SPECIES_CODE %in% c(68560, 68580, 69322, 69323)) |>
    dplyr::select(HAULJOIN, SPECIES_CODE, SEX, SIZE, SHELL_CONDITION, CALCULATED_WEIGHT, WEIGHT, SAMPLING_FACTOR) |>
    dplyr::mutate(FREQUENCY = 1) |>
    dplyr::inner_join(
      dplyr::select(
        all_hauls, HAULJOIN, VESSEL, CRUISE, HAUL, MATCHUP
      ),
      by = "HAULJOIN"
    )

  crab_fish <- dplyr::bind_rows(crab, fish_lengths)
  
  # Identify haul pairs to use for selectivity analysis based on minimum sample size ----
  selectivity_flag <- dplyr::select(crab_fish, HAULJOIN, MATCHUP, SPECIES_CODE, FREQUENCY) |>
    dplyr::group_by(HAULJOIN, MATCHUP, SPECIES_CODE) |>
    dplyr::summarise(N_MEASURED = sum(FREQUENCY), .groups = "keep") |>
    dplyr::inner_join(sratio::species_code_label(x = "all") |> # Set minimum sample size for selectivity analysis
                        dplyr::select(SPECIES_CODE, MIN_SAMPLE_SIZE),
                      by = "SPECIES_CODE") |>
    dplyr::mutate(USE_FOR_SELECTIVITY = N_MEASURED >= MIN_SAMPLE_SIZE) |>
    dplyr::select(-MIN_SAMPLE_SIZE) |>
    dplyr::inner_join(all_hauls[c("HAULJOIN", "MATCHUP", "TREATMENT")],
                      by = c("HAULJOIN", "MATCHUP")) |>
    dplyr::ungroup() |>
    dplyr::group_by(MATCHUP, SPECIES_CODE) |>
    dplyr::summarise(USE_FOR_SELECTIVITY = sum(USE_FOR_SELECTIVITY), .groups = "keep") |>
    dplyr::mutate(USE_FOR_SELECTIVITY = USE_FOR_SELECTIVITY > 1)
  
  catch <- dplyr::inner_join(catch, 
                             selectivity_flag, 
                             by = c("SPECIES_CODE", "MATCHUP"))
  
  crab_fish <- dplyr::inner_join(crab_fish, 
                                 selectivity_flag, 
                                 by = c("SPECIES_CODE", "MATCHUP")) |>
    dplyr::select(HAULJOIN, 
                  MATCHUP, 
                  USE_FOR_SELECTIVITY,
                  VESSEL, 
                  CRUISE, 
                  HAUL, 
                  SPECIES_CODE, 
                  SEX, 
                  SIZE,
                  WEIGHT,
                  CALCULATED_WEIGHT,
                  SHELL_CONDITION,
                  SAMPLING_FACTOR, 
                  FREQUENCY)
  
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

  save(data_1530, file = here::here("data", "data_1530.rda"), compress = "xz")
  
  crab_fish |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::filter(USE_FOR_SELECTIVITY) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY), .groups = "keep") |>
    dplyr::inner_join(sratio::species_code_label(x = "all"),
                      by = "SPECIES_CODE") |>
    dplyr::arrange(YEAR) |>
    write.csv(file = here::here("analysis", "15_30", "plots", "sample_sizes_1530.csv"), row.names = FALSE)
  
  crab_fish |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::filter(USE_FOR_SELECTIVITY) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY), .groups = "keep") |>
    dplyr::inner_join(sratio::species_code_label(x = "all"),
                      by = "SPECIES_CODE") |>
    dplyr::ungroup() |>
    dplyr::arrange(YEAR) |>
    tidyr::pivot_wider(names_from = "YEAR", values_from = "n", values_fill = 0) |>
    dplyr::mutate(COMMON_NAME = sratio::species_code_label(SPECIES_CODE, type = "common_name")) |>
    write.csv(file = here::here("analysis", "15_30", "plots", "sample_sizes_wide_1530.csv"), row.names = FALSE)
  
  
  crab_fish |>
    dplyr::mutate(YEAR = floor(CRUISE/100)) |>
    dplyr::group_by(SPECIES_CODE, YEAR) |>
    dplyr::summarise(n = sum(FREQUENCY), .groups = "keep") |>
    dplyr::inner_join(sratio::species_code_label(x = "all"),
                      by = "SPECIES_CODE") |>
    dplyr::arrange(YEAR) |>
    write.csv(file = here::here("analysis", "15_30", "plots", "sample_sizes_no_filter_1530.csv"), row.names = FALSE)

  cat("---- Retrieving data from gapindex ----\n")
  
  gapindex_1530 <- vector(mode = "list", length = length(species_codes))
  
  setNames(object = gapindex_1530, nm = species_codes)
  
  for(ii in 1:length(species_codes)) {
    
    # Retrieve survey data for selected species
    gapindex_1530[[ii]] <- 
      gapindex::get_data(year_set = 1982:2024, 
                         survey_set = "EBS", 
                         spp_codes = species_codes[ii],
                         pull_lengths = TRUE, 
                         channel = channel)
    
  }
  
  save(gapindex_1530, file = here::here("data", "gapindex_1530.rda"), compress = "xz")
  
}

get_data(species_codes = c(21740, 21720, 10210, 10261, 10110, 10112, 10115, 10120, 10130, 10140, 10180, 10200, 10285, 471, 68560, 68580, 69322))



