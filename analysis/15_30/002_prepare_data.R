library(sratio)

# Treatment levels
treatments <- factor(c(15,30))

# Load built-in data sets
catch_df <- sratio::data_1530$catch |>
  dplyr::filter(USE_FOR_SELECTIVITY)

haul_df <- sratio::data_1530$haul

# Change species codes to include sexes (needed for analyses - could do a lot of refactoring instead)
crab_by_sex <- sratio::data_1530$size |>
  dplyr::filter(SPECIES_CODE %in% c(68560, 68580, 69322)) |>
  dplyr::mutate(SPECIES_CODE = as.numeric(paste0(SPECIES_CODE, SEX)))

size_df <- dplyr::bind_rows(sratio::data_1530$size, crab_by_sex) |>
  dplyr::filter(USE_FOR_SELECTIVITY)

# Data setup ---------------------------------------------------------------------------------------
dat <- dplyr::inner_join(
  haul_df,
  size_df,
  by = c("CRUISE", "MATCHUP", "HAULJOIN", "VESSEL", "HAUL")) |>
  dplyr::mutate(SIZE = dplyr::if_else(!is.na(LENGTH), LENGTH, WIDTH)) |>
  dplyr::select(HAULJOIN, CRUISE, MATCHUP, TREATMENT, AREA_SWEPT_KM2, SPECIES_CODE, SIZE, FREQUENCY, SAMPLING_FACTOR) |> 
  dplyr::mutate(SIZE = round(SIZE),
                FREQUENCY = FREQUENCY * SAMPLING_FACTOR) |>
  dplyr::group_by(HAULJOIN, CRUISE, MATCHUP, SPECIES_CODE, TREATMENT, AREA_SWEPT_KM2, SIZE) |>
  dplyr::summarise(FREQUENCY = sum(FREQUENCY),
                   .groups = "keep") |>
  dplyr::mutate(CPUE_N_KM2 = FREQUENCY/AREA_SWEPT_KM2) |>
  dplyr::group_by(SPECIES_CODE, TREATMENT, SIZE) |>
  dplyr::summarise(FREQUENCY = sum(FREQUENCY),
                   CPUE = sum(CPUE_N_KM2))

dat_by_cruise <- dplyr::inner_join(
  haul_df,
  size_df,
  by = c("CRUISE", "MATCHUP", "HAULJOIN", "VESSEL", "HAUL")) |>
  dplyr::mutate(SIZE = dplyr::if_else(!is.na(LENGTH), LENGTH, WIDTH)) |>
  dplyr::select(HAULJOIN, CRUISE, MATCHUP, TREATMENT, AREA_SWEPT_KM2, SPECIES_CODE, SIZE, FREQUENCY, SAMPLING_FACTOR) |> 
  dplyr::mutate(SIZE = round(SIZE),
                FREQUENCY = FREQUENCY * SAMPLING_FACTOR) |>
  dplyr::group_by(HAULJOIN, CRUISE, MATCHUP, SPECIES_CODE, TREATMENT, AREA_SWEPT_KM2, SIZE) |>
  dplyr::summarise(FREQUENCY = sum(FREQUENCY),
                   .groups = "keep") |>
  dplyr::mutate(CPUE_N_KM2 = FREQUENCY/AREA_SWEPT_KM2) |>
  dplyr::group_by(SPECIES_CODE, TREATMENT, SIZE, CRUISE) |>
  dplyr::summarise(FREQUENCY = sum(FREQUENCY),
                   CPUE = sum(CPUE_N_KM2))

ggplot() +
  geom_bar(data = dat,
                 mapping = aes(x = SIZE, 
                               y = FREQUENCY, 
                               fill = factor(TREATMENT),
                               group = SIZE),
           stat = "identity") +
  scale_fill_colorblind(name = "Treatment") +
  facet_wrap(~SPECIES_CODE, scales = "free") +
  theme_bw()

ggplot() +
  geom_bar(data = test,
           mapping = aes(x = SIZE, 
                         y = CPUE, 
                         fill = factor(TREATMENT),
                         group = SIZE),
           stat = "identity") +
  scale_fill_colorblind(name = "Treatment") +
  facet_wrap(~SPECIES_CODE, scales = "free") +
  theme_bw()


species_codes <- unique(dat$SPECIES_CODE)

sel_dat <- dplyr::filter(dat, SPECIES_CODE == 21720) |>
  dplyr::select(-CPUE) |>
  tidyr::pivot_wider(id_cols = c("SPECIES_CODE", "SIZE"),
                     values_from = "FREQUENCY", 
                     values_fill = 0,
                     names_from = "TREATMENT",
                     names_prefix = "FREQUENCY_") |>
  dplyr::mutate(FREQUENCY_15 = round(FREQUENCY_15),
                FREQUENCY_30 = round(FREQUENCY_30)) |>
  dplyr::mutate(TOTAL_FREQUENCY = FREQUENCY_15 + FREQUENCY_30)

test <- mgcv::gam(FREQUENCY_15/TOTAL_FREQUENCY ~ s(SIZE, bs = "tp"), 
          weights = TOTAL_FREQUENCY, 
          family = binomial(link = "logit"), 
          data = sel_dat)

plot(test)

fit_df <- data.frame(SIZE = min(sel_dat$SIZE):max(sel_dat$SIZE))
fit_df$fit <- predict(test, newdata = fit_df, type = "response")

ggplot() +
  geom_path(data = fit_df, 
             mapping = aes(x = SIZE, y = fit))


sel_dat <- dplyr::filter(dat, SPECIES_CODE == 21720) |>
  dplyr::select(-FREQUENCY) |>
  tidyr::pivot_wider(id_cols = c("SPECIES_CODE", "SIZE"),
                     values_from = "CPUE", 
                     values_fill = 0,
                     names_from = "TREATMENT",
                     names_prefix = "CPUE_") |>
  dplyr::mutate(CPUE_15 = round(CPUE_15),
                CPUE_30 = round(CPUE_30)) |>
  dplyr::mutate(TOTAL_CPUE = CPUE_15 + CPUE_30)

test <- mgcv::gam(CPUE_15/TOTAL_CPUE ~ s(SIZE, bs = "tp"), 
                  weights = TOTAL_CPUE, 
                  family = binomial(link = "logit"), 
                  data = sel_dat)

plot(test)

fit_df <- data.frame(SIZE = min(sel_dat$SIZE):max(sel_dat$SIZE))
fit_df$fit <- predict(test, newdata = fit_df, type = "response")

ggplot() +
  geom_path(data = fit_df, 
            mapping = aes(x = SIZE, y = fit)) +
  geom_hline(yintercept = 0.5) +
  scale_y_continuous(limits = c(0,1))
