library(ggplot2)
library(dplyr)

# Load data ----

cpue_other <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_other.rds"))
cpue_1998 <- readRDS(here::here("analysis", "somerton_2002", "data", "cpue_1998.rds"))


haul_summary_table <- 
  dplyr::bind_rows(
    cpue_other,
    cpue_1998
  ) |>
  dplyr::group_by(
    common_name,
    SEX,
    YEAR
  ) |>
  dplyr::summarise(
    n = n(),
    mean_cpue_15 = format(round(mean(CPUE_NO_KM2_15)), big.mark = ","),
    min_cpue_15 = format(round(min(CPUE_NO_KM2_15)), big.mark = ","),
    max_cpue_15 = format(round(max(CPUE_NO_KM2_15)), big.mark = ","),
    mean_cpue_30 = format(round(mean(CPUE_NO_KM2_30)), big.mark = ","),
    min_cpue_30 = format(round(min(CPUE_NO_KM2_30)), big.mark = ","),
    max_cpue_30 = format(round(max(CPUE_NO_KM2_30)), big.mark = ",")
  ) |>
  dplyr::mutate(
    `Species/sex` = paste0(common_name, " (", SEX, ")"),
    label = paste0(
      n, 
      "\r", mean_cpue_15, " (", min_cpue_15, "-", max_cpue_15, ")",
      "\r", mean_cpue_30, " (", min_cpue_30, "-", max_cpue_30, ")"
      ),
  ) |>
  dplyr::ungroup() |>
  dplyr::select(`Species/sex`, YEAR, label) |>
  tidyr::pivot_wider(
    values_from = label,
    names_from = YEAR,
    values_fill = "0\r-\r-"
  )


write.csv(
  haul_summary_table,
  file = here::here("analysis", "somerton_2002", "plots", "cpue_summary_table.csv"),
  row.names = FALSE
  )
