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

crab_observed_biomass <- vector(mode = "list", length = nrow(crab_species_surveys))
crab_simulated_biomass <- vector(mode = "list", length = nrow(crab_species_surveys))

for(ii in 1:nrow(crab_species_surveys)) {
  
  cat(paste0(ii, " ", crab_species_surveys$species[ii], " ", crab_species_surveys$region[ii]," ", Sys.time(), "\n"))
  
  # Load a list where each object in the list contains stations to drop
  station_draws <- readRDS(
    file = here::here("analysis", "effort_reduction", "output", 
                      sub_dir, paste0(crab_species_surveys$region[ii], "_removed_stations.rds"))
  )
  
  # Set years
  survey_years <- if(crab_species_surveys$region[ii] == "EBS") {1987:2024} else {2010:2024}
  
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
      region = crab_species_surveys$region[ii],
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
        species = crab_species_surveys$species[ii],
        region = crab_species_surveys$region[ii],
        crab_category = "all_categories"
      ) |>
      dplyr::mutate(iter = jj)
  }
  
  crab_simulated_biomass[[ii]] <- do.call(rbind, simulated_reduced)
  
}

crab_biomass_results <- do.call(rbind, crab_simulated_biomass)
crab_observed_results <- do.call(rbind, crab_observed_biomass)

saveRDS(
  object = crab_biomass_results, 
  file = here::here("analysis", 
                    "effort_reduction", 
                    "output", 
                    sub_dir, 
                    paste0(
                      paste(
                        unique(crab_species_surveys$region), 
                        collapse = "_"), 
                      "_crab_biomass_district.rds")
  ),
  compress = "xz"
)

saveRDS(
  object = crab_observed_results, 
  file = here::here("analysis", 
                    "effort_reduction", 
                    "output", 
                    sub_dir, 
                    paste0(
                      paste(
                        unique(crab_species_surveys$region), 
                        collapse = "_"
                      ), "_crab_biomass_observed.rds")
  ),
  compress = "xz"
)
