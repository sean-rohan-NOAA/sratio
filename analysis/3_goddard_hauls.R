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

get_goddard <- function(){
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
                  HAUL < 267) |>
    dplyr::filter(
      (VESSEL == 88 & HAUL %in% c(192:196, 204, 205, 209, 211, 212, 215, 217, 223, 225, 228, 229, 234, 235, 237, 239, 245:246, 258:261))|
        (VESSEL == 89 & HAUL %in% c(207:211, 217, 218, 222, 224, 225, 228, 230:232, 234, 236, 239, 240, 246, 247, 251:252, 262, 263))) |> #removed V89 H 253, 260; V88 H 247, 254 (bc out of station grid)
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
  saveRDS(hauls_1995, here::here("data", "hauls_1995.rds"))
  
}

#get_goddard()

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
