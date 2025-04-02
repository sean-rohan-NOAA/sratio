# Plot effort reduction scenarios
library(sratio)
library(akgfmaps)
library(ggridges)
library(cowplot)

survey_opts <- c("EBS", "NBS") # EBS or NBS
survey_set <- survey_opts[1]

plot_area_ids <- if(survey_set == "EBS") {c(99900, 99901)} else {99902}
loop_area_names <- if(survey_set == "EBS") {c("EBS Standard Plus NW", "EBS Standard")} else{"NBS"}
survey_years <- if(survey_set == "EBS") {1987:2024} else {2010:2024}
survey_definition_id <- ifelse(survey_set == "EBS", 98, 143)

station_loss_pct <- c(15, 25, 35, 50)

# Load species list
akfin_species <- 
  read.csv(
    file = here::here("analysis", "effort_reduction", "data", "akfin_biomass_taxa_by_survey.csv")
  ) |>
  dplyr::filter(SURVEY_DEFINITION_ID == survey_definition_id)

# Setup area ID look-up table for plotting
area_id_lookup <- 
  data.frame(AREA_ID = plot_area_ids,
             AREA_NAME = loop_area_names)

# Combine area ID and species data frames
species_area_lookup <- 
  dplyr::cross_join(
    akfin_species, 
    area_id_lookup
    )

# Load observed fish values
fish_observed <- 
  readRDS(
    file = 
      here::here(
        "analysis", 
        "effort_reduction", 
        "output", 
        station_loss_pct[1], 
        paste0(survey_set, "_fish_observed.rds")
      )
  ) 

fish_observed$biomass_subarea <- 
  fish_observed$biomass_subarea |>
  dplyr::inner_join(species_area_lookup, by = c("SURVEY_DEFINITION_ID", "SPECIES_CODE", "AREA_ID"))

# Load observed crab values

crab_observed <- 
  readRDS(
    file = 
      here::here(
        "analysis", 
        "effort_reduction",
        "output",
        station_loss_pct[1],
        paste0(survey_set, "_crab_biomass_observed.rds"))
  ) |>
  dplyr::select(
    SPECIES, 
    YEAR, 
    REGION, 
    DISTRICT, 
    CATEGORY, 
    FULL_ABUNDANCE = ABUNDANCE,
    FULL_ABUNDANCE_CV = ABUNDANCE_CV, 
    FULL_ABUNDANCE_CI = ABUNDANCE_CI, 
    FULL_BIOMASS_MT = BIOMASS_MT, 
    FULL_BIOMASS_MT_CV = BIOMASS_MT_CV,
    FULL_BIOMASS_MT_CI = BIOMASS_MT_CI, 
    FULL_BIOMASS_LBS = BIOMASS_LBS, 
    FULL_BIOMASS_LBS_CV = BIOMASS_LBS_CV
  )

fish_results <- vector(mode = "list", length = length(station_loss_pct))
crab_results <- vector(mode = "list", length = length(station_loss_pct))

for(ii in 1:length(station_loss_pct)) {
  
  fish_results[[ii]] <- 
    do.call(
      rbind,
      readRDS(
        here::here(
          "analysis", 
          "effort_reduction", 
          "output", 
          station_loss_pct[ii], 
          paste0(survey_set, "_fish_biomass_subarea.rds"))
      )
    )
  
  fish_results[[ii]]$station_loss_pct <- station_loss_pct[ii]
  
  crab_results[[ii]] <- 
      readRDS(
        here::here(
          "analysis", 
          "effort_reduction", 
          "output", 
          station_loss_pct[ii], 
          paste0(survey_set, "_crab_biomass_district.rds"))
    )
  
  crab_results[[ii]]$station_loss_pct <- station_loss_pct[ii]
  
}

fish_results <- do.call(rbind, fish_results)
crab_results <- do.call(rbind, crab_results)

# Calculate fish CVs
fish_results$CV <- sqrt(fish_results$BIOMASS_VAR)/fish_results$BIOMASS_MT

# Setup fish plots and table data
fish_reduced_pct_change <- 
  fish_results |>
  dplyr::inner_join(
    fish_observed$biomass_subarea |>
      dplyr::select(SPECIES_CODE, 
                    AREA_ID,
                    AREA_NAME,
                    COMMON_NAME,
                    YEAR,
                    FULL_BIOMASS_MT = BIOMASS_MT, 
                    FULL_BIOMASS_VAR = BIOMASS_VAR,
                    FULL_CV = CV),
    by = c("SPECIES_CODE", "AREA_ID", "YEAR")
  ) |>
  dplyr::mutate(
    PCT_CV = (CV-FULL_CV)/FULL_CV*100,
    DIFF_CV = CV-FULL_CV,
    PCT_BIOMASS_MT = (BIOMASS_MT-FULL_BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_MT = (BIOMASS_MT-FULL_BIOMASS_MT)/FULL_BIOMASS_MT*100
  ) |>
  dplyr::bind_rows(
    fish_observed$biomass_subarea |>
      dplyr::mutate(
        station_loss_pct = 0,
        PCT_CV = NA,
        DIFF_CV = NA,
        PCT_BIOMASS_MT = NA,
        ABS_PCT_BIOMASS_MT = NA
      )
  )

# Summarize fish changes by species and area
fish_change_conf_int <- 
  fish_reduced_pct_change |>
  dplyr::group_by(SPECIES_CODE, COMMON_NAME, AREA_NAME, station_loss_pct) |>
  dplyr::summarise(
    CV_MEAN = mean(CV, na.rm = TRUE),
    CV_LWR = quantile(CV, 0.025, na.rm = TRUE),
    CV_UPR = quantile(CV, 0.975, na.rm = TRUE),
    PCT_BIOMASS_MT_LWR = quantile(PCT_BIOMASS_MT, 0.025, na.rm = TRUE),
    PCT_BIOMASS_MT_UPR = quantile(PCT_BIOMASS_MT, 0.975, na.rm = TRUE)
  )

# Make fish plots
for(jj in 1:nrow(species_area_lookup)) {
  
  p_title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      paste0(
        species_area_lookup$AREA_NAME[jj],
        " ",
        species_area_lookup$COMMON_NAME[jj]
      ),
      fontface = 'bold', 
      size = 14
    )
  
  p_dat <- 
    dplyr::filter(
      fish_reduced_pct_change, 
      AREA_ID == species_area_lookup$AREA_ID[jj], 
      SPECIES_CODE == species_area_lookup$SPECIES_CODE[jj]
    )
  
  p_biomass <- 
    ggplot() +
    geom_density_ridges(
      data = p_dat,
      mapping = 
        aes(x = PCT_BIOMASS_MT, 
            y = factor(station_loss_pct), 
            fill = station_loss_pct), 
      alpha = 0.6
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_discrete(name = "Station Loss (%)") +
    # scale_x_continuous(name = expression("Biomass change (%), "*(100%*%over(S[R]-S[O], S[O])))) +
    scale_x_continuous(name = "Biomass change (%)") +
    scale_fill_viridis_c() +
    theme_ridges(center_axis_labels = TRUE) +
    theme(legend.position = "none")
  
  p_cv <- 
    ggplot() +
    geom_density_ridges(
      data = p_dat,
      mapping = 
        aes(x = CV, 
            y = factor(station_loss_pct), 
            fill = station_loss_pct), 
      alpha = 0.6, 
      calc_ecdf = TRUE,
      quantiles = 0.5,
      quantile_lines = TRUE
    ) +
    scale_y_discrete(name = "Station Loss (%)") +
    scale_x_continuous(name = "CV") +
    scale_fill_viridis_c() +
    # scale_color_viridis_c() +
    theme_ridges(center_axis_labels = TRUE) +
    theme(legend.position = "none")
  
  p_biomass_cv_grid <-
    cowplot::plot_grid(
      p_title,
      cowplot::plot_grid(
        p_biomass, 
        p_cv,
        align = "hv"
      ),
      rel_heights = c(0.1, 1),
      nrow = 2
    )
  
  png(filename =
        here::here(
          "analysis", 
          "effort_reduction", 
          "plots", 
          paste0(survey_set, "_fish_", 
                 species_area_lookup$AREA_ID[jj], 
                 "_", 
                 species_area_lookup$SPECIES_CODE[jj], 
                 "_diff_all_scenarios.png")
        ),
      width = 6,
      height = 4,
      units = "in",
      res = 300)
  print(p_biomass_cv_grid)
  dev.off()
  
}

# Setup crab plots and table data ----
crab_reduced_pct_change <- 
  dplyr::inner_join(
    crab_observed,
    crab_results,
    by = c(
      "SPECIES", 
      "YEAR", 
      "REGION", 
      "DISTRICT", 
      "CATEGORY"
    )
  ) |>
  dplyr::mutate(
    PCT_ABUNDANCE = (FULL_ABUNDANCE - ABUNDANCE)/FULL_ABUNDANCE*100,
    ABS_PCT_ABUNDANCE = abs((FULL_ABUNDANCE - ABUNDANCE)/FULL_ABUNDANCE)*100,
    PCT_ABUNDANCE_CV = (FULL_ABUNDANCE_CV - ABUNDANCE_CV)/FULL_ABUNDANCE_CV*100,
    DIFF_ABUNDANCE_CV = ABUNDANCE_CV - FULL_ABUNDANCE_CV,
    PCT_BIOMASS_MT = (FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_MT = abs((FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT)*100,
    PCT_BIOMASS_MT_CV = (FULL_BIOMASS_MT_CV - BIOMASS_MT_CV)/FULL_BIOMASS_MT_CV*100,
    DIFF_BIOMASS_MT_CV = BIOMASS_MT_CV - FULL_BIOMASS_MT_CV,
    PCT_BIOMASS_LBS = (FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_LBS = abs((FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT)*100,
    PCT_BIOMASS_LBS_CV = (FULL_BIOMASS_LBS_CV - BIOMASS_LBS_CV)/FULL_BIOMASS_LBS_CV*100,
    DIFF_BIOMASS_LBS_CV = BIOMASS_LBS_CV - FULL_BIOMASS_LBS_CV
  ) |>
  dplyr::bind_rows(
    crab_observed |>
      dplyr::mutate(
        station_loss_pct = 0
      )
  )

# Summarize fish changes by species and area
crab_change_conf_int <- 
  crab_reduced_pct_change |>
  dplyr::group_by(SPECIES, REGION, DISTRICT, CATEGORY, station_loss_pct) |>
  dplyr::summarise(
    CV_MEAN = mean(BIOMASS_MT_CV, na.rm = TRUE),
    CV_LWR = quantile(BIOMASS_MT_CV, 0.025, na.rm = TRUE),
    CV_UPR = quantile(BIOMASS_MT_CV, 0.975, na.rm = TRUE),
    PCT_BIOMASS_MT_LWR = quantile(PCT_BIOMASS_MT, 0.025, na.rm = TRUE),
    PCT_BIOMASS_MT_UPR = quantile(PCT_BIOMASS_MT, 0.975, na.rm = TRUE)
  )

crab_species_area_lookup <- 
  crab_change_conf_int |>
  dplyr::select(SPECIES, REGION, DISTRICT, CATEGORY) |>
  unique()

# Make crab plots
for(kk in 1:nrow(crab_species_area_lookup)) {
  
  p_title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      paste0(
        crab_species_area_lookup$REGION[kk],
        " ",
        crab_species_area_lookup$DISTRICT[kk],
        " ",
        crab_species_area_lookup$SPECIES[kk],
        " ",
        crab_species_area_lookup$CATEGORY[kk]
      ),
      fontface = 'bold', 
      size = 14
    )
  
  p_dat <- 
    crab_reduced_pct_change |>
    dplyr::filter(
      REGION == crab_species_area_lookup$REGION[kk], 
      DISTRICT == crab_species_area_lookup$DISTRICT[kk],
      SPECIES == crab_species_area_lookup$SPECIES[kk],
      CATEGORY == crab_species_area_lookup$CATEGORY[kk]
    )
  
  p_biomass <- 
    ggplot() +
    geom_density_ridges(
      data = p_dat,
      mapping = 
        aes(x = PCT_BIOMASS_MT, 
            y = factor(station_loss_pct), 
            fill = station_loss_pct), 
      alpha = 0.6
    ) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_discrete(name = "Station Loss (%)") +
    # scale_x_continuous(name = expression("Biomass change (%), "*(100%*%over(S[R]-S[O], S[O])))) +
    scale_x_continuous(name = "Biomass change (%)") +
    scale_fill_viridis_c() +
    theme_ridges(center_axis_labels = TRUE) +
    theme(legend.position = "none")
  
  p_cv <- 
    ggplot() +
    geom_density_ridges(
      data = p_dat,
      mapping = 
        aes(x = BIOMASS_MT_CV, 
            y = factor(station_loss_pct), 
            fill = station_loss_pct), 
      alpha = 0.6, 
      calc_ecdf = TRUE,
      quantiles = 0.5,
      quantile_lines = TRUE
    ) +
    scale_y_discrete(name = "Station Loss (%)") +
    scale_x_continuous(name = "CV") +
    scale_fill_viridis_c() +
    # scale_color_viridis_c() +
    theme_ridges(center_axis_labels = TRUE) +
    theme(legend.position = "none")
  
  p_biomass_cv_grid <-
    cowplot::plot_grid(
      p_title,
      cowplot::plot_grid(
        p_biomass, 
        p_cv,
        align = "hv"
      ),
      rel_heights = c(0.1, 1),
      nrow = 2
    )
  
  png(filename =
        here::here(
          "analysis", 
          "effort_reduction", 
          "plots", 
          paste0(survey_set, 
                 "_crab_", 
                 crab_species_area_lookup$SPECIES[kk],
                 " ",
                 crab_species_area_lookup$DISTRICT[kk],
                 " ",
                 crab_species_area_lookup$CATEGORY[kk],
                 "_diff_all_scenarios.png")
        ),
      width = 6,
      height = 4,
      units = "in",
      res = 300)
  print(p_biomass_cv_grid)
  dev.off()
  
}
