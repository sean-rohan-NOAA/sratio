library(sratio)

# Treatment levels
treatments <- factor(c(15,30))

# Load built-in data sets
catch_df <- sratio::data_1530$catch |>
  dplyr::filter(CRUISE %in% use_cruises,
                USE_FOR_SELECTIVITY)

haul_df <- sratio::data_1530$haul |>
  dplyr::filter(CRUISE %in% use_cruises)

size_df <- sratio::data_1530$size |>
  dplyr::filter(CRUISE %in% use_cruises,
                USE_FOR_SELECTIVITY)

# Data setup ---------------------------------------------------------------------------------------
dat <- dplyr::inner_join(haul_df,
             size_df,
             by = c("CRUISE", "MATCHUP", "HAULJOIN", "VESSEL", "HAUL")) |>
  dplyr::select(HAULJOIN, SPECIES_CODE, MATCHUP, CRUISE, FREQUENCY, LENGTH, WIDTH, TREATMENT, DISTANCE_FISHED, NET_WIDTH) |>
  dplyr::mutate(AREA_SWEPT_KM2 = DISTANCE_FISHED * NET_WIDTH / 1000,
                MATCHUP = MATCHUP) |>
  dplyr::group_by(HAULJOIN, SPECIES_CODE, MATCHUP, CRUISE, LENGTH, WIDTH, TREATMENT, AREA_SWEPT_KM2) |>
  dplyr::summarize(FREQUENCY = sum(FREQUENCY), .groups = "keep") |>
  as.data.frame()

# Select measurement to use for size (length or width)
dat$SIZE <- ifelse(!is.na(dat$LENGTH), dat$LENGTH, dat$WIDTH)

# Expand raw length-frequency data to haul-level counts
sampling_factor <- dat |>
  dplyr::group_by(HAULJOIN, SPECIES_CODE) |>
  dplyr::summarize(LEN_COUNT = sum(FREQUENCY), .groups = "keep") |>
  dplyr::inner_join(catch_df, by = c("HAULJOIN", "SPECIES_CODE")) |>
  dplyr::mutate(SAMPLING_FACTOR = NUMBER_FISH/LEN_COUNT) |> # sampling factor
  dplyr::select(HAULJOIN, SPECIES_CODE, SAMPLING_FACTOR, NUMBER_FISH)

dat <- dplyr::inner_join(dat, sampling_factor, by = c("HAULJOIN", "SPECIES_CODE"))

# Expand size-frequency based on sampling factor calculated from total count in catch
dat$FREQ_EXPANDED <- dat$FREQUENCY * dat$SAMPLING_FACTOR

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
    dplyr::summarise(FREQ_EXPANDED = sum(FREQ_EXPANDED), .groups = "keep") |>
    dplyr::bind_rows(dat_binned)
    
  sratio_dat <- sel_dat |>
    dplyr::group_by(HAULJOIN, SIZE_BIN, SPECIES_CODE, MATCHUP) |>
    dplyr::summarise(FREQ_EXPANDED = sum(FREQ_EXPANDED), .groups = "keep") |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = SIZE_BIN, values_from = FREQ_EXPANDED, values_fill = 0)
  
  sratio_dat <- tidyr::pivot_longer(sratio_dat, 
                                    cols = 4:ncol(sratio_dat), 
                                    names_to = "SIZE_BIN", 
                                    values_to = "FREQ_EXPANDED") |>
    dplyr::inner_join(dplyr::select(haul_df, HAULJOIN, TREATMENT, AREA_SWEPT_KM2), by = "HAULJOIN")

  sratio_dat <- sratio_dat |> 
    dplyr::mutate(TREATMENT_COL = paste0("N_", TREATMENT)) |>
    dplyr::select(-HAULJOIN, -TREATMENT, -AREA_SWEPT_KM2) |> 
    tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = FREQ_EXPANDED) |>
    dplyr::mutate(SIZE_BIN = as.numeric(SIZE_BIN)) |>
    dplyr::arrange(MATCHUP, SIZE_BIN) |> 
    dplyr::inner_join(sratio_dat |> 
                        dplyr::mutate(TREATMENT_COL = paste0("AREA_SWEPT_KM2_", TREATMENT)) |>
                        dplyr::select(MATCHUP, TREATMENT_COL, AREA_SWEPT_KM2) |>
                        unique() |>
                        tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = AREA_SWEPT_KM2),
                      by = "MATCHUP"
                      
    )
  
  dat_sratio <- dplyr::bind_rows(dat_sratio, sratio_dat)
  
}

saveRDS(object = dat_binned, file = here::here("output", "catch_at_length_1530.rds"))
saveRDS(object = dat_sratio, file = here::here("output", "n_by_treatment_1530.rds"))
