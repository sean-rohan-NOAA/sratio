library(sratio)

species_code = 10110

bootstrap_results_path <- list.files(here::here("output"), 
                                     recursive = TRUE, 
                                     pattern = "sratio_bootstrap_results_", 
                                     full.names = TRUE)

bootstrap_results <- readRDS(bootstrap_results_path[grep(pattern = species_code, bootstrap_results_path)]) |>
  dplyr::group_by(SIZE_BIN) |>
  dplyr::summarise(s12 = median(s12))

size_30 <- readRDS(here::here("output", "catch_at_length_1530.rds")) |>
  dplyr::filter(TREATMENT == 30, SPECIES_CODE == species_code) |>
  dplyr::ungroup()

size_15 <- readRDS(here::here("output", "catch_at_length_1530.rds")) |>
  dplyr::filter(TREATMENT == 15, SPECIES_CODE == species_code) |>
  dplyr::ungroup()

size_15 <- dplyr::select(size_30, MATCHUP, AREA_SWEPT_KM_30 = AREA_SWEPT_KM2) |> 
  unique() |>
  dplyr::inner_join(size_15) |>
  dplyr::inner_join(bootstrap_results)

size_15$PREDICTED_FREQUENCY <- size_15$FREQ_EXPANDED * size_15$AREA_SWEPT_KM_30/size_15$AREA_SWEPT_KM2 * size_15$s12

size_30 <- dplyr::select(size_15, MATCHUP, PREDICTED_FREQUENCY, SIZE_BIN, SPECIES_CODE) |>
  dplyr::full_join(size_30)

ggplot(data = size_30,
       mapping = aes(x = SIZE_BIN,  y = PREDICTED_FREQUENCY-FREQ_EXPANDED)) +
  geom_boxplot(mapping = aes(group = SIZE_BIN)) +
  geom_smooth()


effort1
effort2
size1


