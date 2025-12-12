library(sratio)

gear_codes <- 
  data.frame(
    GEAR = c("PNE", "83-112"),
    TREATMENT = factor(c(172, 44))
  )

species_codes <-
  data.frame(
    SPECIES_CODE = c(21740, 21720, 20510, 30060, 10115, 10110, 10112, 10130, 471, 68580),
    COMMON_NAME = c("walleye pollock", "Pacific cod", "sablefish", "Pacific ocean perch", "Greenland turbot", "arrowtooth flounder", "Kamchatka flounder", "flathead sole", "Alaska skate", "snow crab")
  )


mean(data_ss$haul$BOTTOM_DEPTH[data_ss$haul$BOTTOM_DEPTH > 0])
range(data_ss$haul$BOTTOM_DEPTH[data_ss$haul$BOTTOM_DEPTH > 0])


fit_mgcv <- function(x, species_code) {
  
  x_sel <- dplyr::filter(x, SPECIES_CODE == species_code)
  
  mod <- mgcv::gam(log(CPUE_KGKM2_172+1) ~ s(log(CPUE_KGKM2_44+1), bs = "tp"), family = gaussian(), data = x)

  newdata <- 
    data.frame(
      SPECIES_CODE = species_code,
      CPUE_KGKM2_44 = seq(min(log(x_sel$CPUE_KGKM2_44+1)), max(log(x_sel$CPUE_KGKM2_44+1)), 0.05)
    )
  
  fit <- predict(mod, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$fit <- fit$fit
  newdata$se.fit <- fit$se.fit
  
  return(
    list(model = mod,
         fit = newdata)
  )
  
}


cpue_data <- 
  dplyr::inner_join(
    data_ss$catch,
    dplyr::select(
      data_ss$haul, VESSEL, CRUISE, HAUL, MATCHUP, TREATMENT, AREA_SWEPT_KM2
    )
  ) |>
  dplyr::inner_join(gear_codes) |>
  dplyr::mutate(
    CPUE_KGKM2 = WEIGHT/AREA_SWEPT_KM2,
    CPUE_NOKM2 = NUMBER_FISH/AREA_SWEPT_KM2,
    COMMON_NAME = factor(SPECIES_CODE, levels = species_codes$SPECIES_CODE, labels = species_codes$COMMON_NAME)) |>
  dplyr::select(-NUMBER_FISH, -HAULJOIN, -VESSEL, -HAUL, -AREA_SWEPT_KM2, -CRUISE, -WEIGHT)


cpue_wide <- cpue_data |>
  dplyr::select(-GEAR) |>
  tidyr::pivot_wider(
    values_from = c("CPUE_KGKM2", "CPUE_NOKM2"),  
    names_from = TREATMENT, 
    values_fill = 0
  )


cpue_data |>
  dplyr::group_by(COMMON_NAME, GEAR) |>
dplyr::summarise(
 n_positive = sum(CPUE_NOKM2 > 0)
) |> 
  dplyr::mutate(
    p_positive = n_positive/40 * 100,
    text = format(p_positive, n_small = 1)
  ) |>
  dplyr::select(-n_positive, -p_positive) |>
  tidyr::pivot_wider(names_from = GEAR, values_from = text)


ggplot() +
  geom_smooth(
    data = cpue_wide,
    mapping = aes(x = log(CPUE_KGKM2_44+1), y = log(CPUE_KGKM2_172+1))
  ) +
  geom_point(
    data = cpue_wide,
    mapping = aes(x = log(CPUE_KGKM2_44+1), y = log(CPUE_KGKM2_172+1))
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~COMMON_NAME, scales = "free") +
  scale_x_continuous(name = expression(ln(CPUE[83-112]+1)*' '*(kg%.%km^-2))) +
  scale_y_continuous(name = expression(ln(CPUE[PNE]+1)*' '*(kg%.%km^-2))) +
  theme_bw()


ggplot() +
  geom_point(
    data = cpue_wide,
    mapping = aes(x = log(CPUE_NOKM2_44+1), y = log(CPUE_NOKM2_172+1))
  ) +
  geom_smooth(
    data = cpue_wide,
    mapping = aes(x = log(CPUE_NOKM2_44+1), y = log(CPUE_NOKM2_172+1))
  )+
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~SPECIES_CODE, scales = "free") +
  theme_bw()



dplyr::filter(cpue_data, SPECIES_CODE == 10112, `172` == 0)
dplyr::filter(cpue_data, SPECIES_CODE == 10112, `44` == 0)
dplyr::filter(cpue_data, SPECIES_CODE == 10115, `44` == 0)