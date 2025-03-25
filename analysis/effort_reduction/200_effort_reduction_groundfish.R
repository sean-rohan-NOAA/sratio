# Effect of losing stations on index data products
# 
# This analysis examines the impact of excluding a proportion of survey stations in the EBS
# on design-based indices of abundance for groundfish and crab. Stations are excluded proportionally 
# from each groundfish subarea. Design-based index products are then calculated for each target 
# taxa.
#
# Created by Sean Rohan
# March 21, 2025

library(sratio) # 2025.03.21
library(akgfmaps) # 4.0.3
library(gapindex) # 3.0.2
library(crabpack) # 1.0.0

# Control pars
prop_drop <- 0.25
n_iter <- 5
seed <- 1337
survey_set <- "EBS" # EBS or NBS
survey_years <- if(survey_set == "EBS") {1987:2024} else {2010:2024}
survey_definition_id <- ifelse(survey_set == "EBS", 98, 143)

# Setup output directories
sub_dir <- gsub(".*\\.", "", prop_drop)
dir.create(path = here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set), 
           showWarnings = FALSE, 
           recursive = TRUE)
dir.create(path = here::here("analysis", "effort_reduction", "output", sub_dir, survey_set), 
           showWarnings = FALSE, 
           recursive = TRUE)

# Load species list
akfin_species <- 
  read.csv(
    file = here::here("analysis", "effort_reduction", "data", "akfin_biomass_taxa_by_survey.csv")
    ) |>
  dplyr::filter(SURVEY_DEFINITION_ID == survey_definition_id)

# Function to remove stations randomly from the EBS/NBS grid in the event of effort reduction
effred_remove_stations <- 
  function(gapdata, stations) {
    
    remove_hauljoin <- gapdata$haul$HAULJOIN[gapdata$haul$STATIONID %in% stations]
    
    gapdata$haul <- gapdata$haul[!(gapdata$haul$HAULJOIN %in% remove_hauljoin), ]
    gapdata$catch <- gapdata$catch[!(gapdata$catch$HAULJOIN %in% remove_hauljoin), ]
    gapdata$specimen <- gapdata$specimen[!(gapdata$specimen$HAULJOIN %in% remove_hauljoin), ]
    
    if(!is.null(gapdata$size)) {
      gapdata$size <- gapdata$size[!(gapdata$size$HAULJOIN %in% remove_hauljoin), ]
    }
    
    return(gapdata)
    
  }

# Connect 
channel <- gapindex::get_connected()

# Assign stations to subareas ----
region_layers <- akgfmaps::get_base_layers(select.region = tolower(survey_set))

# for(jj in 1:nrow(akfin_species)) {

# Calculate baseline data products ----
dat <- 
  gapindex::get_data(
    year_set = survey_years,
    survey_set = survey_set,
    spp_codes = akfin_species$SPECIES_CODE,
    pull_lengths = TRUE,
    channel = channel
  )

observed <- make_gapindex(gapdata = dat)

station_subareas <- 
  sf::st_centroid(region_layers$survey.grid) |>
  sf::st_intersection(region_layers$survey.strata) |>
  sf::st_drop_geometry() |>
  dplyr::select(STATION, STRATUM) |>
  dplyr::inner_join(
    dplyr::filter(
      dat$stratum_groups, 
      AREA_ID %in% 1:9,
      DESIGN_YEAR == max(dat$stratum_groups$DESIGN_YEAR)), 
    by = "STRATUM")

observed$biomass_subarea$CV <- sqrt(observed$biomass_subarea$BIOMASS_VAR)/observed$biomass_subarea$BIOMASS_MT

# Draw removed proportion of stations randomly by subarea and generate products

station_draws <- vector(mode = "list", length = n_iter)
biomass_subarea_results <- vector(mode = "list", length = n_iter)
biomass_stratum_results <- vector(mode = "list", length = n_iter)
agecomp_results <- vector(mode = "list", length = n_iter)
sizecomp_results <- vector(mode = "list", length = n_iter)


start_time <- Sys.time()
set.seed(seed)

for(ii in 1:n_iter) {
  
  print(ii)
  
  # Draw stations
  station_draws[[ii]] <- 
    (station_subareas |>
       dplyr::group_by(AREA_ID) |>
       dplyr::slice_sample(prop = prop_drop))$STATION
  
  # Remove stations
  sel_gapdata <- 
    effred_remove_stations(
      gapdata = dat,
      stations = station_draws[[ii]]
    )
  
  # Generate products
  gapindex_results <-
    make_gapindex(gapdata = sel_gapdata) 
  
  # Add outputs to lists
  biomass_subarea_results[[ii]] <- gapindex_results[['biomass_subarea']] |>
    dplyr::mutate(iter = ii)
  
  biomass_stratum_results[[ii]] <- gapindex_results[['biomass_stratum']] |>
    dplyr::mutate(iter = ii)
  
  if(!is.null(gapindex_results$alk)) {
    sizecomp_results[[ii]] <- gapindex_results[['sizecomp_subarea']] |>
      dplyr::mutate(iter = ii)
    
    agecomp_results[[ii]] <- gapindex_results[['agecomp_region']] |>
      dplyr::mutate(iter = ii)
  }
  
}

# Combine outputs into a data frame
reduced_sampling_results <- do.call(rbind, biomass_subarea_results)

# Calculate CVs
reduced_sampling_results$CV <- sqrt(reduced_sampling_results$BIOMASS_VAR)/reduced_sampling_results$BIOMASS_MT
end_time <- Sys.time()

difftime(end_time, start_time)

# }
# Save groundfish outputs
saveRDS(
  object = observed, 
  file = here::here("analysis", "effort_reduction", "output", sub_dir, paste0(survey_set, "_fish_observed.rds")),
  compress = "xz"
)

saveRDS(
  object = station_draws, 
  file = here::here("analysis", "effort_reduction", "output", sub_dir, paste0(survey_set, "_removed_stations.rds")),
  compress = "xz"
)

saveRDS(
  object = biomass_subarea_results, 
  file = here::here("analysis", "effort_reduction", "output", sub_dir, paste0(survey_set, "_fish_biomass_subarea.rds")),
  compress = "xz"
)

saveRDS(
  object = biomass_stratum_results, 
  file = here::here("analysis", "effort_reduction", "output", sub_dir, paste0(survey_set, "_fish_biomass_stratum.rds")),
  compress = "xz"
)

saveRDS(
  object = agecomp_results, 
  file = here::here("analysis", "effort_reduction", "output", sub_dir, paste0(survey_set, "_fish_agecomp.rds")),
  compress = "xz"
)

saveRDS(
  object = sizecomp_results, 
  file = here::here("analysis", "effort_reduction", "output", sub_dir, paste0(survey_set, "_fish_sizecomp.rds")),
  compress = "xz"
)

# Calculate summary statistics -----

reduced_pct_change <- 
  reduced_sampling_results |>
  dplyr::inner_join(
    ebs_observed$biomass_subarea |>
      dplyr::select(SPECIES_CODE, 
                    AREA_ID,
                    YEAR,
                    FULL_BIOMASS_MT = BIOMASS_MT, 
                    FULL_BIOMASS_VAR = BIOMASS_VAR,
                    FULL_CV = CV),
    by = c("YEAR", "SPECIES_CODE", "AREA_ID")
  ) |>
  dplyr::mutate(
    PCT_CV = (CV-FULL_CV)/FULL_CV*100,
    DIFF_CV = CV-FULL_CV,
    PCT_BIOMASS_MT = (BIOMASS_MT-FULL_BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_MT = (BIOMASS_MT-FULL_BIOMASS_MT)/FULL_BIOMASS_MT*100
  ) |>
  dplyr::filter(AREA_ID %in% c(99900, 99901))


reduced_pct_change_year <- 
  reduced_pct_change |>
  dplyr::group_by(SPECIES_CODE, AREA_ID, YEAR) |>
  dplyr::summarise(
    MEAN_PCT_CV = mean(PCT_CV, na.rm = TRUE),
    MEAN_DIFF_CV = mean(DIFF_CV, na.rm = TRUE),
    MEAN_ABS_PCT_BIOMASS_MT = mean(abs(PCT_BIOMASS_MT), na.rm = TRUE),
    MEAN_PCT_BIOMASS_MT = mean(PCT_BIOMASS_MT, na.rm = TRUE)
  )

reduced_pct_change_overall <- 
  reduced_pct_change |>
  dplyr::group_by(SPECIES_CODE, AREA_ID) |>
  dplyr::summarise(
    MEAN_PCT_CV = mean(PCT_CV, na.rm = TRUE),
    MEAN_DIFF_CV = mean(DIFF_CV, na.rm = TRUE),
    MEAN_ABS_PCT_BIOMASS_MT = mean(abs(PCT_BIOMASS_MT), na.rm = TRUE),
    MEAN_PCT_BIOMASS_MT = mean(PCT_BIOMASS_MT, na.rm = TRUE)
  )


# Loop through Standard and Standard Plus NW 
loop_area_id <- c(99900, 99901)
loop_area_names <- c("EBS Standard Plus NW", "EBS Standard")


for(kk in 1:length(loop_area_id)) {
  
  # Plot: Biomass timeseries ----
  
  p_biomass <- 
    ggplot() +
    geom_line(data = dplyr::filter(
      reduced_sampling_results, 
      AREA_ID == loop_area_id[kk]),
      mapping = aes(
        x = YEAR, 
        y = BIOMASS_MT, 
        color = paste0("Drop ", prop_drop*100, "%"), 
        group = factor(iter)),
      alpha = 0.3) +
    geom_line(data = 
                dplyr::filter(
                  ebs_observed$biomass_subarea, 
                  AREA_ID == loop_area_id[kk]),
              mapping = 
                aes(
                  x = YEAR, 
                  y = BIOMASS_MT, 
                  color = "Observed"),
              linewidth = 1.02
    ) +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Biomass (mt)") +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE, scales = "free")
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0("p_", survey_set, "_", loop_area_id[kk], "_biomass_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_biomass)
  dev.off()
  
  # Plot: CV timeseries ----
  p_cv <- 
    ggplot() +
    geom_line(data = dplyr::filter(
      reduced_sampling_results, 
      AREA_ID == loop_area_id[kk]),
      mapping = aes(
        x = YEAR, 
        y = CV, 
        color = paste0("Drop ", prop_drop*100, "%"), 
        group = factor(iter)),
      alpha = 0.3) +
    geom_line(data = 
                dplyr::filter(
                  ebs_observed$biomass_subarea, 
                  AREA_ID == loop_area_id[kk]),
              mapping = 
                aes(
                  x = YEAR, 
                  y = CV, 
                  color = "Observed"),
              linewidth = 1.02
    ) +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "CV", limits = c(0, NA)) +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE)
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0("p_", survey_set, "_", loop_area_id[kk], "_cv_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_cv)
  dev.off()
  
  # Plot: Biomass percent change timeseries ----
  p_biomass_pct <-
    ggplot() +
    geom_jitter(
      data =
        dplyr::filter(
          reduced_pct_change,
          AREA_ID == loop_area_id[kk]
        ),
      mapping = aes(x = YEAR, y = PCT_BIOMASS_MT),
      size = rel(0.2),
      color = "grey50",
      alpha = 0.3, 
      width = 0.25,
      height = 0) +
    geom_line(data = 
                dplyr::filter(
                  reduced_pct_change_year,
                  AREA_ID == loop_area_id[kk]
                ),
              mapping = aes(
                x = YEAR, 
                y = MEAN_PCT_BIOMASS_MT),
              size = 1.05
              ) +
    geom_hline(data = 
                 dplyr::filter(
                   reduced_pct_change_overall,
                   AREA_ID == loop_area_id[kk]
                 ),
               mapping = 
                 aes(yintercept = MEAN_PCT_BIOMASS_MT),
               linetype = 2) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Biomass difference (%)") +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE)
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0("p_", survey_set, "_",  loop_area_id[kk], "_biomass_pct_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_biomass_pct)
  dev.off()
  
  # Plot: Biomass absolute percent change timeseries ----
  p_biomass_abs_pct <-
    ggplot() +
    geom_jitter(
      data =
        dplyr::filter(
          reduced_pct_change,
          AREA_ID == loop_area_id[kk]
        ),
      mapping = aes(x = YEAR, y = ABS_PCT_BIOMASS_MT),
      size = rel(0.3),
      color = "grey50",
      alpha = 0.2, 
      width = 0.25,
      height = 0) +
    geom_line(data = 
                dplyr::filter(
                  reduced_pct_change_year,
                  AREA_ID == loop_area_id[kk]
                ),
              mapping = aes(x = YEAR, y = MEAN_ABS_PCT_BIOMASS_MT)) +
    geom_hline(data = 
                 dplyr::filter(
                   reduced_pct_change_overall,
                   AREA_ID == loop_area_id[kk]
                 ),
               mapping = 
                 aes(yintercept = MEAN_ABS_PCT_BIOMASS_MT),
               linetype = 2) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Mean biomass change (%)", limits = c(0, NA)) +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE)
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0("p_", survey_set, "_",  loop_area_id[kk], "_biomass_abs_pct_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_biomass_abs_pct)
  dev.off()
  
  # Plot: CV change ----
  
  p_cv_change <- 
    ggplot() +
    geom_jitter(
      data =
        dplyr::filter(
          reduced_pct_change,
          AREA_ID == loop_area_id[kk]
        ),
      mapping = aes(x = YEAR, y = DIFF_CV),
      size = rel(0.2),
      color = "grey50",
      alpha = 0.5) +
    geom_line(
      data = 
        dplyr::filter(
          reduced_pct_change_year,
          AREA_ID == loop_area_id[kk]
        ),
      mapping = 
        aes(x = YEAR, 
            y = MEAN_DIFF_CV)
    ) +
    geom_hline(
      data = 
        dplyr::filter(
          reduced_pct_change_overall,
          AREA_ID == loop_area_id[kk]
        ),
      mapping = 
        aes(yintercept = MEAN_DIFF_CV),
      linetype = 2) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = expression('Mean '*CV[reduced]-CV[full])) +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE)
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0("p_", survey_set, "_", loop_area_id[kk], "_cv_diff.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_cv_change)
  dev.off()
  
}
