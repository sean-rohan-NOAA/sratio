library(sratio)

trawl_geometry <-
sratio::data_1530$haul |>
  dplyr::filter(NET_MEASURED == 'Y') |>
  dplyr::group_by(MATCHUP) |>
  dplyr::summarise(n = n()) |>
  dplyr::filter(n >= 2) |>
  dplyr::select(MATCHUP) |>
  dplyr::inner_join(sratio::data_1530$haul) |>
  dplyr::mutate(SCOPE_TO_DEPTH = WIRE_LENGTH/BOTTOM_DEPTH,
                UNIQUE_NET = if_else(is.na(NET_NUMBER), NA, interaction(CRUISE, VESSEL, NET_NUMBER))
                ) |>
  dplyr::mutate(UNIQUE_NET = as.factor(as.character(UNIQUE_NET))) |>
  as.data.frame()

# Spread and height GAMs ----
# 104 hauls where we could track individual nets

width_gam_0 <- mgcv::gam(NET_WIDTH ~ + s(UNIQUE_NET, bs = "re"), 
                         data = trawl_geometry)

width_gam_1 <- mgcv::gam(NET_WIDTH ~ s(SCOPE_TO_DEPTH, bs = "tp") + s(UNIQUE_NET, bs = "re"), 
                         data = trawl_geometry)

width_gam_2 <- mgcv::gam(NET_WIDTH ~ s(SCOPE_TO_DEPTH, bs = "tp") + s(TOW_SPEED_KNOTS, bs = "tp") + s(UNIQUE_NET, bs = "re"), 
                         data = trawl_geometry)

width_gam_3 <- mgcv::gam(NET_WIDTH ~ TREATMENT + s(SCOPE_TO_DEPTH, bs = "tp") + s(UNIQUE_NET, bs = "re"), 
                         data = trawl_geometry)

height_gam_0 <- mgcv::gam(NET_HEIGHT ~ + s(UNIQUE_NET, bs = "re"), 
                          data = trawl_geometry)

height_gam_1 <- mgcv::gam(NET_HEIGHT ~ s(SCOPE_TO_DEPTH, bs = "tp") + s(UNIQUE_NET, bs = "re"), 
                          data = trawl_geometry)

height_gam_2 <- mgcv::gam(NET_HEIGHT ~ s(SCOPE_TO_DEPTH, bs = "tp") + s(TOW_SPEED_KNOTS, bs = "tp") + s(UNIQUE_NET, bs = "re"), 
                          data = trawl_geometry)

height_gam_3 <- mgcv::gam(NET_HEIGHT ~ TREATMENT + s(SCOPE_TO_DEPTH, bs = "tp") + s(UNIQUE_NET, bs = "re"), 
                          data = trawl_geometry)


make_aic_table <- function(x) {
  
  output <- data.frame()
  
  for(ii in 1:length(x)) {
    
    output <- rbind(output, 
                    data.frame(formula = paste0(as.character(x[[ii]]$formula)[2], "~", as.character(x[[ii]]$formula)[3]),
               edf = round(sum(x[[ii]]$edf) + x[[ii]]$nsdf, 2),
               AIC = round(x[[ii]]$aic, 2),
               deltaAIC = NA,
               deviance_explained = round((x[[ii]]$null.deviance-x[[ii]]$deviance)/x[[ii]]$null.deviance * 100, 1))
    )
    
    output$deltaAIC <- round(output$AIC - min(output$AIC), 2)
    
  }
  
  return(output)
  
}


dplyr::bind_rows(
  data.frame(formula = "Spread"),
  make_aic_table(x = list(width_gam_0, width_gam_1, width_gam_2, width_gam_3)),
  data.frame(formula = "Height"),
  make_aic_table(x = list(height_gam_0, height_gam_1, height_gam_2, height_gam_3)),
  data.frame(formula = "Bottom contact")
) |>
  write.csv(file = here::here("analysis", "15_30", "plots", "spread_height_gam_table.csv"), row.names = FALSE)


table_trawl_geometry <- sratio::data_1530$haul |>
  dplyr::group_by(TREATMENT) |>
  dplyr::summarise(N_HAULS = n(),
                   MEAN_TOW_SPEED_KNOTS = format(mean(TOW_SPEED_KNOTS), nsmall = 2, digits = 2),
                   SD_TOW_SPEED_KNOTS = format(sd(TOW_SPEED_KNOTS), nsmall = 2, digits = 2),
                   MIN_TOW_SPEED_KNOTS = format(min(TOW_SPEED_KNOTS), nsmall = 2, digits = 2),
                   MAX_TOW_SPEED_KNOTS = format(max(TOW_SPEED_KNOTS), nsmall = 2, digits = 2),
                   MEAN_NET_HEIGHT = format(mean(NET_HEIGHT), nsmall = 2, digits = 2),
                   SD_NET_HEIGHT = format(sd(NET_HEIGHT), nsmall = 2, digits = 2),
                   MIN_NET_HEIGHT = format(min(NET_HEIGHT), nsmall = 2, digits = 2),
                   MAX_NET_HEIGHT = format(max(NET_HEIGHT), nsmall = 2, digits = 2),
                   MEAN_NET_WIDTH = format(mean(NET_WIDTH), nsmall = 2, digits = 2),
                   SD_NET_WIDTH = format(sd(NET_WIDTH), nsmall = 2, digits = 2),
                   MIN_NET_WIDTH = format(min(NET_WIDTH), nsmall = 2, digits = 2),
                   MAX_NET_WIDTH = format(max(NET_WIDTH), nsmall = 2, digits = 2)
                   )

data.frame(Treatment = table_trawl_geometry$TREATMENT,
           Hauls = table_trawl_geometry$N_HAULS,
           Spread = paste0(table_trawl_geometry$MEAN_NET_WIDTH, 
                           " (", table_trawl_geometry$MIN_NET_WIDTH, 
                           "-", table_trawl_geometry$MAX_NET_WIDTH, ")"),
           Height = paste0(table_trawl_geometry$MEAN_NET_HEIGHT, 
                           " (", 
                           table_trawl_geometry$MIN_NET_HEIGHT, 
                           "-", table_trawl_geometry$MAX_NET_HEIGHT, ")"),
           Speed = paste0(table_trawl_geometry$MEAN_TOW_SPEED_KNOTS, 
                          " (", table_trawl_geometry$MIN_TOW_SPEED_KNOTS, 
                          "-", table_trawl_geometry$MAX_TOW_SPEED_KNOTS, ")")
           ) |>
write.csv(file = here::here("analysis", "15_30", "plots", "speed_spread_height.csv"), row.names = FALSE)

