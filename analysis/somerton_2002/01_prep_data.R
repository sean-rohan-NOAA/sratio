library(nortest)
library(xlsx)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(dplyr)
library(here)

# Somerton's archived data

somerton_catch <-
  rbind(
    read.table(file = here::here("analysis", "somerton_2002", "data", "Cb2.txt"),
               header = TRUE, na.strings = ".") |>
      dplyr::mutate(SPECIES_CODE = 68560),
    read.table(file = here::here("analysis", "somerton_2002", "data", "Co2.txt"),
               header = TRUE, na.strings = ".") |>
      dplyr::mutate(SPECIES_CODE = 68580),
    read.table(file = here::here("analysis", "somerton_2002", "data", "RK2.txt"),
               header = TRUE, na.strings = ".") |>
      dplyr::mutate(SPECIES_CODE = 69322)
  ) |>
  dplyr::mutate(CATCH15F = CATCH15T-CATCH15M,
                CATCH30F = CATCH30T-CATCH30M) |>
  dplyr::select(TOW_PAIR = OBS, SPECIES_CODE, CATCH15M, CATCH30M, CATCH15F, CATCH30F) |>
  tidyr::pivot_longer(cols = c("CATCH15M", "CATCH30M", "CATCH15F", "CATCH30F")) |>
  dplyr::mutate(TREATMENT = factor(ifelse(stringr::str_detect(name, "30"), 30, 15)),
                SEX = factor(ifelse(stringr::str_detect(name, "F"), "F", "M"))) |>
  dplyr::select(-name) |>
  tidyr::pivot_wider(values_from = "value", names_from = "TREATMENT", names_prefix = "COUNT_")

somerton_effort <-
  read.table(file = here::here("analysis", "somerton_2002", "data", "Cb1.txt"),
             header = TRUE, na.strings = ".") |>
  dplyr::select(TOW_PAIR = OBS, EFFORT15, EFFORT30, VESSEL) |>
  tidyr::pivot_longer(cols = c("EFFORT15", "EFFORT30")) |>
  dplyr::mutate(TREATMENT = factor(ifelse(stringr::str_detect(name, "30"), 30, 15))) |>
  dplyr::select(-name) |>
  dplyr::mutate(value = value/100,
                VESSEL = factor(VESSEL)) |>
  tidyr::pivot_wider(values_from = "value", names_from = "TREATMENT", names_prefix = "AREA_SWEPT_KM2_")

cpue_1998 <-
  dplyr::inner_join(somerton_catch, somerton_effort) |>
  dplyr::mutate(
    YEAR = 1998,
    CPUE_NO_KM2_15 = COUNT_15 / AREA_SWEPT_KM2_15,
    CPUE_NO_KM2_30 = COUNT_30 / AREA_SWEPT_KM2_30,
    LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
    LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
    CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
    CPUE_RATIO = CPUE_NO_KM2_15/CPUE_NO_KM2_30,
    COMBINED_COUNT = ceiling(COUNT_30 + COUNT_15),
    PROP_15 = CPUE_NO_KM2_15/(CPUE_NO_KM2_15+CPUE_NO_KM2_30),
    EFFORT_RATIO = AREA_SWEPT_KM2_15/AREA_SWEPT_KM2_30,
    COUNT_15 = ceiling(COUNT_15),
    COUNT_30 = ceiling(COUNT_30),
    common_name = sratio::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0)

saveRDS(object = cpue_1998, file = here::here("analysis", "somerton_2002", "data", "cpue_1998.rds"))


# 1995 + 2021-2024 experiments

cpue_other <-
  sratio::data_1530$size |>
  dplyr::filter(SPECIES_CODE %in% c(68580, 68560, 69322)) |>
  dplyr::mutate(
    SEX = ifelse(SEX == 1, "M", "F"),
    MATCHUP = paste0("O", MATCHUP)
  ) |>
  dplyr::group_by(VESSEL, CRUISE, HAUL, SPECIES_CODE, SEX) |>
  dplyr::summarise(COUNT = ceiling(sum(SAMPLING_FACTOR))) |>
  dplyr::inner_join(
    dplyr::select(
      sratio::data_1530$haul,
      VESSEL, CRUISE, HAUL, AREA_SWEPT_KM2, TREATMENT, YEAR, TOW_PAIR = MATCHUP
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(YEAR != 1998) |>
  dplyr::select(-VESSEL, -CRUISE, -HAUL) |>
  tidyr::pivot_wider(values_from = c("COUNT", "AREA_SWEPT_KM2"), names_from = TREATMENT, values_fill = 0) |>
  dplyr::mutate(
    CPUE_NO_KM2_15 = COUNT_15 / AREA_SWEPT_KM2_15,
    CPUE_NO_KM2_30 = COUNT_30 / AREA_SWEPT_KM2_30,
    LOG_CPUE_NO_KM2_30 = log(CPUE_NO_KM2_30),
    LOG_CPUE_NO_KM2_15 = log(CPUE_NO_KM2_15),
    CPUE_LOG_RATIO = log(CPUE_NO_KM2_15/CPUE_NO_KM2_30),
    CPUE_RATIO = CPUE_NO_KM2_15/CPUE_NO_KM2_30,
    COMBINED_COUNT = ceiling(COUNT_30 + COUNT_15),
    PROP_15 = CPUE_NO_KM2_15/(CPUE_NO_KM2_15+CPUE_NO_KM2_30),
    EFFORT_RATIO = AREA_SWEPT_KM2_15/AREA_SWEPT_KM2_30,
    common_name = sratio::species_code_label(SPECIES_CODE, type = "common_name")) |>
  dplyr::filter(CPUE_NO_KM2_15 > 0, CPUE_NO_KM2_30 > 0)

saveRDS(object = cpue_other, file = here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))