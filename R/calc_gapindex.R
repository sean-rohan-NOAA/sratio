#' Function to calculate design-based products from gapindex tables and age-length keys
#' 
#' Wrapper around gapindex functions to calcualte survey biomass, abundance, age comps, and size comps.
#' 
#' @param racebase_tables data object created from gapindex::get_data()
#' @param alk age-length key (dataframe) created from gapindex::calc_alk()
#' @param subset_fields Should all of the fields be returned for biomass, agecomp, and sizecomp data frames?
#' @param fill_NA_method passed to gapindex::calc_sizecomp_stratum()
#' @export
#' @import gapindex

calc_gapindex <- function(racebase_tables, alk = NULL, subset_fields = FALSE, fill_NA_method = "AIGOA") {
  
  if(is.null(alk)) {
    alk <- gapindex::calc_alk(
      racebase_tables = racebase_tables
    )
  }
  
  cpue <- 
    gapindex::calc_cpue(
      racebase_tables = racebase_tables
    )
  
  biomass_stratum <- 
    suppressWarnings(
      gapindex::calc_biomass_stratum(
        racebase_tables = racebase_tables,
        cpue = cpue
      )
    )
  
  biomass_subarea <- 
    suppressWarnings(
      gapindex::calc_biomass_subarea(
        racebase_tables = racebase_tables,
        biomass_strata = biomass_stratum)
    )
  
  sizecomp_stratum <- 
    gapindex::calc_sizecomp_stratum(
      racebase_tables = racebase_tables,
      racebase_cpue = cpue, 
      racebase_stratum_popn = biomass_stratum,
      spatial_level = "stratum",
      fill_NA_method = fill_NA_method
    )
  
  sizecomp_subarea <- 
    gapindex::calc_sizecomp_subarea(
      racebase_tables = racebase_tables,
      size_comps = sizecomp_stratum
    )
  
  agecomp_stratum <- 
    gapindex::calc_agecomp_stratum(
      racebase_tables = racebase_tables, 
      alk = alk, 
      size_comp = sizecomp_stratum
    )
  
  agecomp_region <- 
    suppressWarnings(
      gapindex::calc_agecomp_region(racebase_tables = racebase_tables,
                                    age_comps_stratum = agecomp_stratum
      )
    )
  
  if(subset_fields) {
    return(
      list(biomass = 
             dplyr::filter(biomass_subarea, AREA_ID %in% c(99900, 99901)) |>
             dplyr::select(SPECIES_CODE, AREA_ID, YEAR, BIOMASS_MT, BIOMASS_VAR, POPULATION_COUNT, POPULATION_VAR),
           sizecomp = 
             dplyr::filter(sizecomp_subarea, AREA_ID %in% c(99900, 99901)) |>
             dplyr::select(SPECIES_CODE, AREA_ID, YEAR, LENGTH_MM, SEX, POPULATION_COUNT),
           agecomp = dplyr::filter(agecomp_region, AREA_ID %in% c(99900, 99901)) |>
             dplyr::select(SPECIES_CODE, AREA_ID, YEAR, SEX, AGE, LENGTH_MM_MEAN)
      )
    )
  }
  
  return(
    list(biomass = biomass_subarea,
         sizecomp = sizecomp_stratum,
         agecomp = agecomp_region
    )
  )
  
}