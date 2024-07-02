library(sratio)
library(akgfmaps)


channel <- get_connected(schema = "AFSC")

# Retrieve data from the 15/30 hauls from Pam Goddard's thesis
# Goddard, P. D. 1997. The effects of tow duration and subsampling on CPUE, species composition and 
#   length distributions of bottom trawl survey catches. University of Washington. 119 pp.
hauls_1995 <- RODBC::sqlQuery(channel = channel,
                              query = "select * from racebase.haul
                                        where cruise = 199501
                                        and vessel in (88, 89)
                                        and haul_type = 7
                                        and performance >= 0") |>
  dplyr::filter(!is.na(NET_WIDTH),
                !is.na(DISTANCE_FISHED),
                PERFORMANCE >= 0) |>
  dplyr::filter(DISTANCE_FISHED > 0,
                START_TIME > as.POSIXct("1995-07-24 PDT")) |>
  dplyr::filter(!(HAUL %in% c(218) & VESSEL == 88),
                HAUL < 267)

hauls_1995_sf <- sf::st_as_sf(hauls_1995, coords = c("START_LONGITUDE", "START_LATITUDE"), crs = "WGS84")


map_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "WGS84", fix.invalid.geom = FALSE)

sf::st_make_valid(map_layers$bathymetry)

ggplot() +
  geom_sf(data = hauls_1995_sf,
          mapping = aes(color = factor(VESSEL)))

joined_hauls <- hauls_1995_sf |>
  dplyr::rename(OLD_STATIONID = STATIONID) |>
  sf::st_join(map_layers$survey.grid)

apply(joined_hauls, MARGIN = 1, FUN = gsub(pattern = ""))

joined_hauls$STATIONID


View(dplyr::select(joined_hauls, HAULJOIN, VESSEL, CRUISE, HAUL, OLD_STATIONID, STATIONID))


ggplot() +
  geom_sf(data = map_layers$survey.grid |>
            dplyr::filter(STATIONID %in% c("AZ0504", "A-06", "A-05", "A-04", "Z-05", "A-07", "B-08", "B-07", "A-06", "B-06")),
          mapping = aes(fill = factor(STATIONID))) +
  geom_sf(data = joined_hauls |>
            dplyr::filter(OLD_STATIONID %in% c("000A-06", "000Z-05", "B-8")))

ggplot() +
  geom_sf(data = map_layers$survey.grid |>
            dplyr::filter(STATIONID %in% c("L-18", "K-18")),
          mapping = aes(fill = factor(STATIONID))) +
  geom_sf(data = joined_hauls |>
            dplyr::filter(OLD_STATIONID %in% c("000L-18", "000K-18")),
          mapping = aes(fill = factor(STATIONID)),
          shape = 21)

ggplot() +
  geom_sf(data = map_layers$survey.grid |>
            dplyr::filter(STATIONID %in% c("B-08", "B-07", "A-06", "B-06")),
          mapping = aes(fill = factor(STATIONID))) +
  geom_sf(data = joined_hauls |>
            dplyr::filter(OLD_STATIONID %in% c("B-8")))

