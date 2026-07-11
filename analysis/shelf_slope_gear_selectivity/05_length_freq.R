# Supplement: Full survey length-frequency versus study length-frequency

library(gapindex)
library(ggthemes)
library(ggridges)

# Load survey size comp data for EBS shelf 2002-2025 and EBS slope 2002-2016
sizecomp_plot_data <- 
  readRDS(
    object = sizecomp_plot_data, 
    file = here::here("analysis", "shelf_slope_gear_selectivity", "sizecomp_plot_data.rds")
  )

# Load catch comparison experiment size comp data
length_comp <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp.rds"))

study_length_comp <- 
  length_comp |>
  dplyr::mutate(
    NUMBER_FISH = round(FREQUENCY*RAISING_FACTOR),
    GEAR = ifelse(GEAR == 44, "83-112", "PNE")
  ) |>
  dplyr::group_by(
    COMMON_NAME, LENGTH, GEAR
  ) |>
  dplyr::summarise(NUMBER_FISH = sum(NUMBER_FISH)) |>
  dplyr::ungroup()

study_length_comp <- 
  study_length_comp |>
  dplyr::group_by(COMMON_NAME, GEAR) |>
  dplyr::summarise(
    TOTAL_COUNT = sum(NUMBER_FISH)
  ) |>
  dplyr::ungroup() |>
  dplyr::inner_join(
    study_length_comp
  ) |>
  dplyr::mutate(PROP = NUMBER_FISH/TOTAL_COUNT) |>
  dplyr::rename(LENGTH_CM = LENGTH,
                SURVEY = GEAR,
                SPECIES_CODE = COMMON_NAME)


(p_survey_vs_study_sizecomp <- 
  ggplot() +
  ggridges::geom_density_ridges(
    data = dplyr::bind_rows(
      sizecomp_plot_data,
      study_length_comp),
    mapping = aes(
      x = LENGTH_CM,
      y = factor(
        SURVEY, levels = c("PNE", "83-112", "EBS Shelf", "EBS Slope")),
      height = PROP,
      fill = SURVEY), 
    stat = "identity",
    color = "grey20",
    alpha = 0.8
  ) +
  facet_wrap(~factor(SPECIES_CODE, levels = taxa_levels), scales = "free") +
  scale_fill_manual(
    name = NULL,
    values = c(
      "PNE" = "#E69F00",
      "83-112" = "#000000",
      "EBS Shelf" = "#56B4E9",
      "EBS Slope" = "#009E73"
    ),
    labels = c(
      "PNE" = "PNE (this study)",
      "83-112" = "83-112 (this study)",
      "EBS Shelf" = "EBS Shelf Survey (2002-2025)",
      "EBS Slope" = "EBS Slope Survey (2002-2016)"
    )
  ) +
  scale_x_continuous(name = "Length (cm)") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9),
        legend.position = "inside",
        legend.position.inside = c(0.70, 0.12),
        axis.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold")))

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "survey_vs_study_sizecomp.png"),
    width = 169, height = 169, units = "mm", res = 300)
print(p_survey_vs_study_sizecomp)
dev.off()
