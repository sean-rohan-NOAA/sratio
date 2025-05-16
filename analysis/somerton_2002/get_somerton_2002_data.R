# Get Somerton data
library(sratio)

crab_species_codes <- c(68560, 68580, 69322)

channel <- sratio::get_connected()

somerton_hauls_1998 <- 
  RODBC::sqlQuery(
    channel = channel,
    query = 
      paste0(
        "SELECT hauljoin, vessel, cruise, haul, duration, distance_fished, net_width 
      FROM racebase.haul
      WHERE
        cruise = 199801
        AND hauljoin IN (", paste(unique(sratio::crab_size_1995_1998$HAULJOIN), collapse = ", "),
        ")"
      ) 
  ) |>
  dplyr::mutate(
    AREA_SWEPT_KM2 = NET_WIDTH/1000*DISTANCE_FISHED,
    TREATMENT = 
      factor(
        plyr::round_any(DURATION, 0.25), 
        levels = c(0.25, 0.5), 
        labels = c(15, 30)
      )
  ) |>
  dplyr::inner_join(
    dplyr::select(
      sratio::otto_key_1998, # Use tow pairs from Somerton et al. (2002)
      VESSEL, CRUISE, HAUL, TOW_PAIR, TOW_BLOCK
    ),
    by = c("VESSEL", "CRUISE", "HAUL")
  ) |>
  dplyr::arrange(TOW_BLOCK)

somerton_crab <- 
  somerton_hauls_1998 |>
  dplyr::select(HAULJOIN, VESSEL, CRUISE, HAUL, AREA_SWEPT_KM2, DURATION, TREATMENT, TOW_PAIR, TOW_BLOCK) |>
  dplyr::inner_join(
    dplyr::filter(sratio::crab_size_1995_1998),
    by = c("HAULJOIN", "VESSEL", "CRUISE", "HAUL")
  ) |>
  dplyr::filter(SPECIES_CODE %in% crab_species_codes,  # RKC, snow, Tanner from 1998
                CRUISE == 199801) |>
  dplyr::mutate(SPECIES_CODE_SEX  = as.numeric(paste0(SPECIES_CODE, SEX)),  # Combine sex and species code
                VESSEL = factor(VESSEL),
                SEX = factor(SEX)) # Make vessel a factor for models

saveRDS(somerton_crab, here::here("analysis", "somerton_2002", "data", "somerton_crab.rds"))
saveRDS(somerton_hauls_1998, here::here("analysis", "somerton_2002", "data", "somerton_hauls_1998.rds"))