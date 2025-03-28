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
sub_dir <- gsub(".*\\.", "", prop_drop)
n_iter <- 100
seed <- 1337
survey_set <- "NBS" # EBS or NBS
survey_years <- if(survey_set == "EBS") {1987:2024} else {2010:2024}
survey_definition_id <- ifelse(survey_set == "EBS", 98, 143)

# Setup output directories
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

# Identify subareas (EBS) or strata (NBS) for each station ----
region_layers <- akgfmaps::get_base_layers(select.region = tolower(survey_set))

if(survey_set == "EBS") {
  
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
  
} 

if(survey_set == "NBS") {
  
  station_subareas <- 
    sf::st_centroid(region_layers$survey.grid) |>
    sf::st_intersection(region_layers$survey.strata) |>
    sf::st_drop_geometry() |>
    dplyr::select(STATION, AREA_ID = STRATUM)
  
}

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
  
  cat(paste0(ii, " ", Sys.time(), "\n"))
  
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

end_time <- Sys.time()

difftime(end_time, start_time)