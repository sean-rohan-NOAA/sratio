# Effect of losing stations on crab products
# 
# This analysis examines the impact of excluding a proportion of survey stations in the EBS
# on design-based indices of abundance. Stations are excluded proportionally from each groundfish 
# stratum.
#
# Created by Sean Rohan
# March 25, 2025

library(crabpack) # 1.0.0

prop_drop <- 0.25
sub_dir <- gsub(".*\\.", "", prop_drop)

crab_species_surveys <- 
  expand.grid(
    species = c("RKC", "BKC", "TANNER", "SNOW", "HYBRID", "HAIR"),
    region = c("EBS", "NBS")
  )

survey_set <- "EBS" # EBS or NBS
survey_years <- if(survey_set == "EBS") {1987:2024} else {2010:2024}
survey_definition_id <- ifelse(survey_set == "EBS", 98, 143)

# Load a list where each object in the list contains stations to drop
station_draws <- readRDS(
  file = here::here("analysis", "effort_reduction", "output", 
                    sub_dir, paste0(survey_set, "_removed_stations.rds"))
  )

crab_observed_biomass <- vector(mode = "list", length = nrow(crab_species_surveys))
crab_simulated_biomass <- vector(mode = "list", length = nrow(crab_species_surveys))

for(ii in 1:nrow(crab_species_surveys)) {
  
  # Retrieve crab data from 1987-2024
  crab_data <- crabpack::get_specimen_data(
    species = crab_species_surveys$species[ii], 
    region = crab_species_surveys$region[ii],
    years = survey_years,
    channel = "API"
  )
  
  # Calculate observed biomass and abundance
  crab_observed_biomass[[ii]] <- 
    crabpack::calc_bioabund(
      crab_data = crab_data, 
      species = crab_species_surveys$species[ii],
      crab_category = "all_categories"
    )
  
  simulated_reduced <- vector(mode = "list", length = length(station_draws))
  
  for(jj in 1:length(station_draws)) {
    # Drop stations ----
    sel_crab_data <- crab_data
    
    sel_crab_data$haul <- sel_crab_data$haul[!(sel_crab_data$haul$STATION_ID %in% station_draws[[jj]]), ]
    
    sel_crab_data$specimen <- sel_crab_data$specimen[!(sel_crab_data$specimen$STATION_ID %in% station_draws[[jj]]), ]
    
    # Calculate biomass and abundance after dropping stations
    simulated_reduced[[jj]] <- 
      crabpack::calc_bioabund(
        crab_data = sel_crab_data, 
        species = crab_species[ii],
        crab_category = "all_categories"
      ) |>
      dplyr::mutate(iter = jj)
  }
  
  crab_simulated_biomass <- simulated_reduced
  
}



