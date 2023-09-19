# The Goddard experiment didn't truly conduct side-by-side tows. There are clusters of tows and vessels 
# didnt necessarily start towing at the same time. In some cases tow paths may have slightly overlapped.
# Another potential issue is vessel avoidance behaviors because vessels were operating in the same
# location for multiple hauls.
# 
# Should these data be included in the analysis? If so, how should the matchups be done?


library(sratio)
library(akgfmaps)
library(plotly)

map_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

channel <- get_connected(schema = "AFSC")

# Retrieve data from the 15/30 hauls from Pam Goddard's thesis
# Goddard, P. D. 1997. The effects of tow duration and subsampling on CPUE, species composition and 
#   length distributions of bottom trawl survey catches. University of Washington. 119 pp.
hauls_1995 <- RODBC::sqlQuery(channel = channel,
                              query = "select * from racebase.haul
                                        where cruise = 199501
                                        and vessel in (88, 89)
                                        and haul_type = 7") |>
  dplyr::filter(!is.na(NET_WIDTH),
                !is.na(DISTANCE_FISHED),
                PERFORMANCE >= 0) |>
  dplyr::filter(DISTANCE_FISHED > 0,
                START_TIME > as.POSIXct("1995-07-24 PDT")) |>
  dplyr::filter(!(HAUL %in% c(218) & VESSEL == 88),
                HAUL < 267)

# Make an sf object showing haul paths (estimated as a straight line between start and end positions)
hauls_1995_sf <- hauls_1995 |>
  dplyr::select(START_LONGITUDE, END_LONGITUDE, START_LATITUDE, END_LATITUDE, VESSEL, HAUL, STATIONID, START_TIME, DURATION) |>
  tidyr::pivot_longer(cols = c("START_LONGITUDE", "END_LONGITUDE")) |>
  dplyr::rename(LONGITUDE = value) |>
  dplyr::mutate(LATITUDE = ifelse(name == "START_LONGITUDE", START_LATITUDE, END_LATITUDE)) |>
  dplyr::select(-START_LATITUDE, -END_LATITUDE, -name) |>
  sf::st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = "EPSG:4326") |>
  sf::st_transform(crs = "EPSG:3338") |>
  dplyr::group_by(VESSEL, HAUL, STATIONID, START_TIME, DURATION) |>
  dplyr::summarise(n = n()) |>
  dplyr::select(-n) |>
  sf::st_cast(to = "LINESTRING")


# Evaluating which hauls should be removed
plotly::ggplotly(
ggplot() +
  geom_sf(data = hauls_1995_sf,
          mapping = aes(color = factor(VESSEL),
                        text = paste0(HAUL, " ", as.character(START_TIME)))) +
  scale_color_discrete(name = "Vessel")
)

ggplot() +
  geom_sf(data = map_layers$survey.strata,
          fill = NA) +
  geom_sf(data = hauls_1995_sf,
          mapping = aes(color = factor(VESSEL),
                        text = paste0(HAUL, " ", as.character(START_TIME)))) +
  scale_color_discrete(name = "Vessel")
