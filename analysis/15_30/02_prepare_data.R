library(sratio)

# Treatment levels
treatments <- factor(c(15,30))

# Load built-in data sets
catch_df <- sratio::data_1530$catch |>
  dplyr::filter(USE_FOR_SELECTIVITY)

haul_df <- sratio::data_1530$haul

# Change species codes to include sexes (needed for analyses - could do a lot of refactoring instead)
crab_sizes <- sratio::data_1530$size |>
  dplyr::filter(SPECIES_CODE %in% c(68560, 68580, 69322)) |>
  dplyr::mutate(SPECIES_CODE = as.numeric(paste0(SPECIES_CODE, SEX)))

fish_sizes <- sratio::data_1530$size |>
  dplyr::filter(!(SPECIES_CODE %in% c(68560, 68580, 69322)))

size_df <- dplyr::bind_rows(crab_sizes, fish_sizes) |>
  dplyr::filter(USE_FOR_SELECTIVITY)

# Data setup ---------------------------------------------------------------------------------------
dat <- dplyr::inner_join(
  haul_df,
  size_df,
  by = c("CRUISE", "MATCHUP", "HAULJOIN", "VESSEL", "HAUL")) |>
  dplyr::mutate(SIZE = dplyr::if_else(!is.na(LENGTH), LENGTH, WIDTH)) |>
  dplyr::select(HAULJOIN, CRUISE, MATCHUP, TREATMENT, AREA_SWEPT_KM2, SPECIES_CODE, SIZE, FREQUENCY, SAMPLING_FACTOR)

species_codes <- unique(dat$SPECIES_CODE)

dat_binned <- data.frame()
dat_sratio <- data.frame()

for(ii in 1:length(species_codes)) {
  
  sel_dat <- dplyr::filter(dat, SPECIES_CODE == species_codes[ii])
  
  sel_dat$SIZE_BIN <- sratio:::make_size_bins(sel_dat$SIZE, species_code = species_codes[ii])
  
  # Only use hauls where at least one individual was caught in both 15 and 30 minute tows
  check_complete <- sel_dat |>
    dplyr::group_by(MATCHUP, TREATMENT) |>
    dplyr::summarize(n = n(), .groups = "keep") |>
    tidyr::pivot_wider(id_cols = MATCHUP, 
                       values_from = n, 
                       names_from = TREATMENT, 
                       values_fill = 0)
  
  check_complete <- check_complete[as.vector(check_complete[, 2] > 0 & check_complete[, 3] > 0), ] 
  
  sel_dat <- dplyr::filter(sel_dat, MATCHUP %in% check_complete$MATCHUP)
  
  dat_binned <- sel_dat |>
    dplyr::group_by(HAULJOIN, MATCHUP, SPECIES_CODE, TREATMENT, AREA_SWEPT_KM2, SIZE_BIN) |>
    dplyr::summarise(SAMPLING_FACTOR = sratio::weighted_mean(x = SAMPLING_FACTOR, w = FREQUENCY),
                     FREQUENCY = sum(FREQUENCY),
                     .groups = "keep") |>
    dplyr::bind_rows(dat_binned)
    
  sratio_dat <- sel_dat |>
    dplyr::group_by(HAULJOIN, MATCHUP, SPECIES_CODE, TREATMENT, AREA_SWEPT_KM2, SIZE_BIN) |>
    dplyr::summarise(FREQUENCY = sum(FREQUENCY),
                     .groups = "keep") |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = SIZE_BIN, values_from = FREQUENCY, values_fill = 0)
  
  sratio_dat <- tidyr::pivot_longer(sratio_dat, 
                                    cols = 6:ncol(sratio_dat), 
                                    names_to = "SIZE_BIN", 
                                    values_to = "FREQUENCY") |>
    dplyr::mutate(SIZE_BIN = as.numeric(SIZE_BIN))
  
  sratio_length <- sratio_dat |> 
    dplyr::mutate(TREATMENT_COL = paste0("N_", TREATMENT)) |>
    dplyr::select(-HAULJOIN, -TREATMENT, -AREA_SWEPT_KM2) |> 
    tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = FREQUENCY)
  
  sratio_area_swept <- sratio_dat |> 
    dplyr::mutate(TREATMENT_COL = paste0("AREA_SWEPT_KM2_", TREATMENT)) |>
    dplyr::select(MATCHUP, TREATMENT_COL, AREA_SWEPT_KM2) |>
    unique() |>
    tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = AREA_SWEPT_KM2)
  
  sratio_sampling_factor <- sel_dat |>
    dplyr::mutate(TREATMENT_COL = paste0("SAMPLING_FACTOR_", TREATMENT)) |>
    dplyr::group_by(MATCHUP, TREATMENT_COL, SPECIES_CODE, SIZE_BIN) |>
    dplyr::summarise(SAMPLING_FACTOR = sratio::weighted_mean(x = SAMPLING_FACTOR, w = FREQUENCY),
                     .groups = "keep") |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = SAMPLING_FACTOR, values_fill = 1)
  
  sratio_dat <- sratio_length |> 
    dplyr::left_join(sratio_area_swept, by = "MATCHUP") |>
    dplyr::left_join(sratio_sampling_factor, by = c("MATCHUP", "SPECIES_CODE", "SIZE_BIN")) |>
    dplyr::mutate(SAMPLING_FACTOR_15 = dplyr::if_else(is.na(SAMPLING_FACTOR_15), 1, SAMPLING_FACTOR_15),
                  SAMPLING_FACTOR_30 = dplyr::if_else(is.na(SAMPLING_FACTOR_30), 1, SAMPLING_FACTOR_30))
  
  dat_sratio <- dplyr::bind_rows(dat_sratio, sratio_dat)
  
}

saveRDS(object = dat_binned, 
        file = here::here("analysis", "15_30", "output", "catch_at_length_1530.rds"))
saveRDS(object = dat_sratio, 
        file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds"))
