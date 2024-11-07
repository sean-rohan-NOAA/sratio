library(gapindex)
library(sratio)

analysis_species_codes <- c(471, 10110, 10130, 10210, 10261, 10285, 21720, 21740)

by_quantile <- TRUE

# Paths to selectivity ratio bootstrap results and gapindex data
fpath_bootstrap_fits <- list.files(path = here::here("analysis", "15_30", "output"),
                                   pattern = "sratio_bootstrap_results",
                                   recursive = TRUE,
                                   full.names = TRUE
)

fpath_gapindex_data <- list.files(path = here::here("analysis", "15_30", "data"),
                                  pattern = "gapindex_data_",
                                  recursive = TRUE,
                                  full.names = TRUE
)

start <- Sys.time()

for(ii in 1:length(analysis_species_codes)) {
  
  # Load gapindex species_data for selected species
  species_data <- readRDS(file = fpath_gapindex_data[grepl(x = fpath_gapindex_data, pattern = analysis_species_codes[ii])])
  
  # Setup size data for simulation
  sim_data <- species_data
  species_size_table <- sim_data$size
  species_size_table$OBSERVED_FREQUENCY <- species_size_table$FREQUENCY
  
  species_catch_table <- sim_data$catch
  species_catch_table$OBSERVED_NUMBER_FISH <- species_catch_table$NUMBER_FISH
  
  # Calculate age-length key
  alk <- gapindex::calc_alk(racebase_tables = species_data)
  
  # Calculate biomass, abundance, length comp, and age comp for the selected species
  index_comps_30_min <- sratio::calc_gapindex(racebase_tables = species_data,
                                             alk = alk,
                                             subset_fields = TRUE)
  
  # Observed sizes
  unique_sizes <- min(sim_data$size$LENGTH/10):max(sim_data$size$LENGTH/10)
  
  # Load bootstrap for selected species
  fits_by_draw <- readRDS(file = fpath_bootstrap_fits[grepl(x = fpath_bootstrap_fits, pattern = analysis_species_codes[ii])]) |>
    dplyr::select(SPECIES_CODE, LENGTH = SIZE_BIN, s12, draw)
  
  fits_by_quantile <- fits_by_draw |>
    dplyr::group_by(LENGTH, SPECIES_CODE) |>
    dplyr::summarise(`0.025` = quantile(s12, 0.025),
                     `0.25` = quantile(s12, 0.25),
                     `0.50` = quantile(s12, 0.5),
                     `0.75` = quantile(s12, 0.75),
                     `0.975` = quantile(s12, 0.975)) |>
    tidyr::pivot_longer(cols = c("0.025", "0.25", "0.50", "0.75", "0.975"),
                        names_to = "quantile",
                        values_to = "s12")
  
  # Sizes encountered during experimental tows
  min_size_fit <- min(fits_by_draw$LENGTH)
  max_size_fit <- max(fits_by_draw$LENGTH)
  
  experiment_sizes <- min_size_fit:max_size_fit
  
  # Sizes that have been observed during the survey that were not encountered during experimental tows
  missing_sizes <- unique_sizes[unique_sizes < min_size_fit | unique_sizes > max_size_fit]
  
  below_min <- missing_sizes[missing_sizes < min_size_fit]
  above_max <- missing_sizes[missing_sizes > max_size_fit] 
  
  # Setup a list to store outputs
  if(by_quantile) {
    unique_levels <- unique(fits_by_quantile$quantile)
  } else {
    unique_levels <- unique(fits_by_draw$draw)
  }

  
  index_comps_15_min <- vector(mode = "list", 
                               length = length(unique_levels))
  
  # Recalculate length composition for each bootstrap sample
  for(jj in 1:length(unique_levels)) {
    
    cat(analysis_species_codes[ii], ": ", jj, "/", length(unique_levels), "\n")
    
    # Select data from draw
    if(by_quantile) {
      sel_fit <-  dplyr::filter(fits_by_quantile, quantile == unique_levels[jj])
      sel_quantile <- unique_levels[jj]
    } else{
      sel_fit <-  dplyr::filter(fits_by_draw, draw == unique_levels[jj]) |>
        dplyr::select(-draw)
      sel_quantile <- NULL
    }

    
    sel_fit <- data.frame(SPECIES_CODE = sel_fit$SPECIES_CODE[ii],
                          LENGTH = c(below_min, above_max),
                          s12 = c(rep(sel_fit$s12[sel_fit$LENGTH == min_size_fit], length(below_min)),  
                                  rep(sel_fit$s12[sel_fit$LENGTH == max_size_fit], length(above_max)))
    ) |>
      dplyr::bind_rows(sel_fit) |>
      dplyr::mutate(LENGTH = LENGTH * 10) |>
      unique() |>
      dplyr::arrange(LENGTH)
    
    sim_data$size <- dplyr::inner_join(species_size_table, 
                                       sel_fit,
                                       by = join_by(SPECIES_CODE, LENGTH))
    
    sim_data$size$FREQUENCY <- sim_data$size$FREQUENCY / sim_data$size$s12
    
    # Estimate total numbers based on size comp
    sim_data$catch <- dplyr::group_by(sim_data$size, HAULJOIN) |>
      dplyr::summarise(FREQUENCY_MULTIPLIER = sum(FREQUENCY)/sum(OBSERVED_FREQUENCY)) |>
      dplyr::ungroup() |>
      dplyr::inner_join(species_catch_table, by = "HAULJOIN") |>
      dplyr::mutate(NUMBER_FISH = NUMBER_FISH * FREQUENCY_MULTIPLIER)
    
    index_comps_15_min[[jj]] <- sratio::calc_gapindex(racebase_tables = sim_data,
                                                      alk = alk,
                                                      subset_fields = TRUE)
    
    index_comps_15_min[[jj]]$biomass$quantile <- as.numeric(sel_quantile)
    index_comps_15_min[[jj]]$sizecomp$quantile <- as.numeric(sel_quantile)
    index_comps_15_min[[jj]]$agecomp$quantile <- as.numeric(sel_quantile)
    
  }
  
  biomass_index_15_min <- do.call(
    rbind, 
    lapply(
      index_comps_15_min, 
      function(x) 
        x$biomass[c("AREA_ID", "YEAR", "BIOMASS_MT", "POPULATION_COUNT", "quantile")]
    )
  )
  
  output <- list(species_data = species_data,
               alk = alk,
               index_comps_30_min = index_comps_30_min, 
               index_comps_15_min = biomass_index_15_min)
  
  saveRDS(object = output, 
          file = 
            here::here("analysis", 
                       "15_30", 
                       "output", 
                       analysis_species_codes[ii], 
                       paste0("sratio_index_comps_15_min_", 
                              analysis_species_codes[ii], 
                              ".rds")
                       )
          )
  print(
  ggplot() +
    geom_point(data = biomass_index_15_min,
               mapping = aes(x = YEAR, y = POPULATION_COUNT, color = "Estimated (15)")) +
    geom_path(data = index_comps_30_min$biomass,
              mapping = aes(x = YEAR, y = POPULATION_COUNT, color = "Observed (30")) +
    ggtitle(label = sratio::species_code_label(x = analysis_species_codes[ii], type = "common_name")) +
    facet_wrap(~AREA_ID)
  )
  
  print(
    ggplot() +
      geom_path(data = biomass_index_15_min[biomass_index_15_min$quantile == 0.5,],
                 mapping = aes(x = YEAR, 
                               y = POPULATION_COUNT,
                               color = "Est. (15)")) +
      geom_path(data = biomass_index_15_min[biomass_index_15_min$quantile == 0.25,],
                mapping = aes(x = YEAR, 
                              y = POPULATION_COUNT,
                              color = "Est. (15)"),
                linetype = 2) +
      geom_path(data = biomass_index_15_min[biomass_index_15_min$quantile == 0.75,],
                mapping = aes(x = YEAR, 
                              y = POPULATION_COUNT,
                              color = "Est. (15)"),
                linetype = 2) +
      geom_ribbon(data = biomass_index_15_min[biomass_index_15_min$quantile %in% 0.75,],
                mapping = aes(x = YEAR, 
                              y = POPULATION_COUNT,
                              color = "Est. (15)"),
                linetype = 2) +
      geom_path(data = index_comps_30_min$biomass,
                mapping = aes(x = YEAR, y = POPULATION_COUNT, color = "Observed (30)")) +
      ggtitle(label = sratio::species_code_label(x = analysis_species_codes[ii], type = "common_name")) +
      facet_wrap(~AREA_ID)
  )
  
  intermediate_time <- Sys.time()
  print(intermediate_time-start)
  
}

end <- Sys.time()
print(end-start)


