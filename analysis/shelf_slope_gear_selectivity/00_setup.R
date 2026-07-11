# Setup data
source(here::here("0_functions.R"))

analysis_species <- 
  data.frame(
    SPECIES_CODE = c(21740, 21720, 10110, 10112, 30060, 10130, 420, 435, 440, 455, 471, 472, 475, 477, 480, 485),
    COMMON_NAME = c(
      "walleye pollock", "Pacific cod", "arrowtooth flounder", "Kamchatka flounder", "Pacific ocean perch", "flathead sole",
      rep("skates", 10))
  )

unique_taxa <- unique(analysis_species$COMMON_NAME)

area_swept <-
  dplyr::select(data_ss$haul, MATCHUP, AREA_SWEPT_KM2, GEAR) |>
  unique() |>
  tidyr::pivot_wider(names_from = "GEAR", values_from = "AREA_SWEPT_KM2", names_prefix = "AREA_SWEPT_KM2_")


area_swept <- dplyr::cross_join(
  dplyr::select(analysis_species, COMMON_NAME) |> unique(),
  area_swept
)

catch_data <- 
  data_ss$catch |>
  dplyr::inner_join(
    analysis_species
  ) |>
  dplyr::group_by(
    HAULJOIN, COMMON_NAME, MATCHUP
  ) |>
  dplyr::summarise(
    WEIGHT = sum(WEIGHT, na.rm = TRUE),
    NUMBER_FISH = sum(NUMBER_FISH, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::inner_join(
    dplyr::select(data_ss$haul, HAULJOIN, GEAR)
  )

catch_data_wide <- catch_data |> 
  dplyr::select(WEIGHT, NUMBER_FISH, COMMON_NAME, MATCHUP, GEAR) |>
  tidyr::pivot_wider(
    names_from = "GEAR", 
    values_from = c("WEIGHT", "NUMBER_FISH"),
    values_fill = 0) |>
  dplyr::inner_join(area_swept)

saveRDS(
  object = catch_data,
  file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "catch_data.rds")
)

saveRDS(
  object = catch_data_wide,
  file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "catch_data_wide.rds")
)

# Unbinned lengths
length_comp <- 
  data_ss$size |>
  dplyr::inner_join(analysis_species) |>
  dplyr::group_by(COMMON_NAME, LENGTH, HAULJOIN, MATCHUP) |>
  dplyr::summarise(
    FREQUENCY = sum(FREQUENCY, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::inner_join(
    dplyr::select(data_ss$haul, HAULJOIN, GEAR, MATCHUP)
  ) |>
  dplyr::select(-HAULJOIN)

# Positive-catch haul pairs
positive_catch_pairs <- 
  length_comp |>
  dplyr::select(COMMON_NAME, MATCHUP, GEAR) |>
  unique() |>
  dplyr::group_by(COMMON_NAME, MATCHUP) |>
  dplyr::summarise(N = n()) |>
  dplyr::filter(N>1) |>
  dplyr::ungroup() |>
  dplyr::select(COMMON_NAME, MATCHUP)

raising_factors <- 
  length_comp |>
  dplyr::group_by(
    COMMON_NAME, GEAR, MATCHUP
  ) |>
  dplyr::summarise(SUM_FREQUENCY = sum(FREQUENCY, na.rm = TRUE)) |> 
  dplyr::ungroup() |>
  dplyr::inner_join(
    catch_data
  ) |>
  dplyr::mutate(RAISING_FACTOR = NUMBER_FISH/SUM_FREQUENCY) |>
  dplyr::select(COMMON_NAME, GEAR, MATCHUP, RAISING_FACTOR)

length_comp <- 
  dplyr::inner_join(
    length_comp, raising_factors
  ) |>
  dplyr::inner_join(
    positive_catch_pairs
  )

length_comp_wide <- length_comp |>
  tidyr::pivot_wider(
    names_from = "GEAR",
    values_from = c("FREQUENCY", "RAISING_FACTOR"),
    values_fill = 0
  ) |>
  dplyr::inner_join(
    area_swept
  )

length_comp_wide$RAISING_FACTOR_172[length_comp_wide$RAISING_FACTOR_172 == 0] <- 1
length_comp_wide$RAISING_FACTOR_44[length_comp_wide$RAISING_FACTOR_44 == 0] <- 1


saveRDS(length_comp, file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp.rds"))

saveRDS(length_comp_wide, file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp_wide.rds"))

# Binned lengths
binned_comp <- 
  data_ss$size |>
  dplyr::mutate(LENGTH_BIN = round_any(LENGTH, accuracy = 4)) |>
  dplyr::inner_join(analysis_species) |>
  dplyr::group_by(COMMON_NAME, LENGTH_BIN, HAULJOIN, MATCHUP) |>
  dplyr::summarise(
    FREQUENCY = sum(FREQUENCY, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::inner_join(
    dplyr::select(data_ss$haul, HAULJOIN, GEAR, MATCHUP)
  ) |>
  dplyr::select(-HAULJOIN)

raising_factors <- 
  binned_comp |>
  dplyr::group_by(
    COMMON_NAME, GEAR, MATCHUP
  ) |>
  dplyr::summarise(SUM_FREQUENCY = sum(FREQUENCY, na.rm = TRUE)) |> 
  dplyr::ungroup() |>
  dplyr::inner_join(
    catch_data
  ) |>
  dplyr::mutate(RAISING_FACTOR = NUMBER_FISH/SUM_FREQUENCY) |>
  dplyr::select(COMMON_NAME, GEAR, MATCHUP, RAISING_FACTOR)

binned_comp <- 
  dplyr::inner_join(
    binned_comp, raising_factors
  ) |>
  dplyr::inner_join(positive_catch_pairs)

binned_comp_wide <- binned_comp |>
  tidyr::pivot_wider(
    names_from = "GEAR",
    values_from = c("FREQUENCY", "RAISING_FACTOR"),
    values_fill = 0
  ) |>
  dplyr::inner_join(
    area_swept
  )

binned_comp_wide$RAISING_FACTOR_172[binned_comp_wide$RAISING_FACTOR_172 == 0] <- 1
binned_comp_wide$RAISING_FACTOR_44[binned_comp_wide$RAISING_FACTOR_44 == 0] <- 1

saveRDS(binned_comp, file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "binned_length_comp.rds"))

saveRDS(binned_comp_wide, file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "binned_length_comp_wide.rds"))

# Bootstrap length samples for each species --- takes a little while

unique_species <- unique(analysis_species$COMMON_NAME)

bootstrap_binned_comp_wide <- vector(mode = "list", length = length(unique_species))
names(bootstrap_binned_comp_wide) <- unique_species

for(ii in 1:length(unique_species)) {
  
  message(Sys.time(), ": ", ii, " ", unique_species[ii])
  
  step1 <- dplyr::filter(
    binned_comp_wide, 
    COMMON_NAME == unique_species[ii])
  
  step2 <- 
    two_stage_bootstrap(
      count1 = step1$FREQUENCY_172[step1$FREQUENCY_172 > 0],
      size1 = step1$LENGTH_BIN[step1$FREQUENCY_172 > 0],
      block1 = step1$MATCHUP[step1$FREQUENCY_172 > 0],
      count2 = step1$FREQUENCY_44[step1$FREQUENCY_44 > 0],
      size2 = step1$LENGTH_BIN[step1$FREQUENCY_44 > 0],
      block2 = step1$MATCHUP[step1$FREQUENCY_44 > 0],
      treatment_name1 = 172,
      treatment_name2 = 44,
      n_draws = 1000, # 1000 bootstrap samples
      seed = 1782850892 # Generated from numeric(Sys.time())
    )  
  
  step3 <- 
    lapply(
      step2,
      function(x, rf, as, common_name) {
        names(x) <- c("GEAR", "LENGTH_BIN", "FREQUENCY", "MATCHUP", "NEW_MATCHUP")
        x$COMMON_NAME <- common_name
        x <- dplyr::left_join(x, rf, by = c("GEAR", "MATCHUP", "COMMON_NAME"))
        x <- tidyr::pivot_wider(x, names_from = "GEAR", values_from = c("FREQUENCY", "RAISING_FACTOR"), values_fill = 0)
        x$RAISING_FACTOR_44[x$RAISING_FACTOR_44 == 0] <- 1
        x$RAISING_FACTOR_172[x$RAISING_FACTOR_172 == 0] <- 1
        x <- dplyr::left_join(x, as, by = c("COMMON_NAME", "MATCHUP"))
        x$ORIGINAL_MATCHUP <- x$MATCHUP
        x$MATCHUP <- x$NEW_MATCHUP
        x <- dplyr::select(x, -NEW_MATCHUP)
        return(x)
      },
      rf = raising_factors,
      as = area_swept,
      common_name = unique_species[ii]
    )
  
  bootstrap_binned_comp_wide[[ii]] <- step3
  
}

saveRDS(
  bootstrap_binned_comp_wide,
  here::here("analysis", "shelf_slope_gear_selectivity", "data", "bootstrap_binned_comp_wide.rds")
)


# Retrieve length-freq from GAP_PRODUCTS (AFSC internal)
library(gapindex) # Github: afsc-gap-products/gapindex

channel <- gapindex::get_connected(check_access = FALSE)

analysis_taxa <- 
  data.frame(
    SPECIES_CODE = c(21740, 21720, 30060, 10110, 10112, 10130, 420, 435, 440, 455, 471, 472, 475, 477, 480, 485),
    GROUP_CODE = c("walleye pollock", "Pacific cod", "Pacific ocean perch", "arrowtooth flounder", "Kamchatka flounder", "flathead sole", rep("skates", 10))
  )

taxa_levels <- unique(analysis_taxa$GROUP_CODE)

bss_dat <- 
  gapindex::get_data(
    year_set = 2002:2025,
    survey_set = c("BSS", "EBS"),
    spp_codes = analysis_taxa,
    pull_lengths = TRUE,
    channel = channel
  )

bss_cpue <- gapindex::calc_cpue(bss_dat)

bss_biomass_stratum <- 
  gapindex::calc_biomass_stratum(
    gapdata = bss_dat, 
    cpue = bss_cpue
  )

bss_biomass_subarea <- 
  gapindex::calc_biomass_subarea(
    gapdata = bss_dat, 
    biomass_stratum = bss_biomass_stratum
  )

bss_sizecomp_stratum <- 
  gapindex::calc_sizecomp_stratum(
    gapdata = bss_dat, 
    cpue = bss_cpue, 
    abundance_stratum = bss_biomass_stratum
  )

bss_sizecomp_subarea <- 
  gapindex::calc_sizecomp_subarea(
    gapdata = bss_dat, 
    sizecomp_stratum = bss_sizecomp_stratum
  )

sizecomp_plot_data <- 
  bss_sizecomp_subarea |>
  dplyr::filter(AREA_ID %in% c(99900, 99905)) |>
  dplyr::group_by(SURVEY_DEFINITION_ID, AREA_ID, SPECIES_CODE, LENGTH_MM) |>
  dplyr::summarize(POPULATION_COUNT = sum(POPULATION_COUNT, na.rm = TRUE)) |>
  dplyr::mutate(LENGTH_CM = LENGTH_MM/10) |>
  dplyr::ungroup()

sizecomp_plot_data <- 
  bss_sizecomp_subarea |>
  dplyr::filter(AREA_ID %in% c(99900, 99905)) |>
  dplyr::group_by(SURVEY_DEFINITION_ID, AREA_ID, SPECIES_CODE) |>
  dplyr::summarize(TOTAL_COUNT = sum(POPULATION_COUNT, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::inner_join(sizecomp_plot_data) |>
  dplyr::mutate(PROP = POPULATION_COUNT/TOTAL_COUNT,
                SURVEY = ifelse(SURVEY_DEFINITION_ID == 78, "EBS Slope", "EBS Shelf"))

saveRDS(object = sizecomp_plot_data, file = here::here("analysis", "shelf_slope_gear_selectivity", "sizecomp_plot_data.rds"))
