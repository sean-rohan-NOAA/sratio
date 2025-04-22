# Prep crab data from Somerton

library(sratio)

channel <- sratio::get_connected()

somerton_crab_1998 <-
  read.csv(here::here("assets", "somerton_crab_spec_1998.csv")) |>
  dplyr::mutate(WIDTH = as.numeric(WIDTH))

crab_spec_1995 <-
  read.csv(here::here("assets", "crab_spec_1995.csv"))

hauls_1998 <-
  RODBC::sqlQuery(
    channel = channel,
    query = "SELECT vessel, cruise, haul, hauljoin, stationid
    FROM racebase.haul
    WHERE cruise = 199801
    AND region = 'BS'
    ")

somerton_crab_1998 <-
  dplyr::left_join(somerton_crab_1998, hauls_1998)

crab_size_1995_1998  <-
  dplyr::bind_rows(
  crab_spec_1995,
  somerton_crab_1998)

save(crab_size_1995_1998, file = here::here("data", "crab_size_1995_1998.rda"))

summary(somerton_crab_1998)

save(crab_size_1995_1998, file = here::here("data", "crab_size_1995_1998.rda"))
