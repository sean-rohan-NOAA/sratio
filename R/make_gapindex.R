#' Generate gapindex results
#'
#' This function is a wrapper function that uses the output of `gapindex::get_data()` to generate selected data product outputs (see below).
#'
#' @param gapdata A data frame containing survey data required for index calculations.
#'
#' @return A list containing the following data tables:
#'   - `cpue`: Catch per unit effort (CPUE) estimates.
#'   - `biomass_stratum`: Biomass estimates at the stratum level.
#'   - `biomass_subarea`: Biomass estimates at the subarea level.
#'   - `alk`: Age-length key (ALK).
#'   - `sizecomp_stratum`: Size composition estimates at the stratum level.
#'   - `sizecomp_subarea`: Size composition estimates at the subarea level.
#'   - `agecomp_stratum`: Age composition estimates at the stratum level.
#'   - `agecomp_region`: Age composition estimates at the regional level.
#'
#' @examples
#' \dontrun{
#' library(gapindex)
#' 
#' gapdata <- 
#' gapindex::get_data(
#'      year_set = 1987:2024, 
#'      survey_set = "EBS", 
#'      spp_codes = 21720, 
#'      pull_lengths = TRUE
#'      )
#' 
#' result <- make_gapindex(gapdata)
#' }
#' @import gapindex
#' @export

make_gapindex <- function(gapdata) {
  
  cpue <- calc_cpue(gapdata = gapdata)
  
  biomass_stratum <- 
    calc_biomass_stratum(
      gapdata = gapdata, 
      cpue = cpue
    )
  
  biomass_subarea <- 
    calc_biomass_subarea(
      gapdata = gapdata, 
      biomass_stratum = biomass_stratum
    )
  
  alk <- calc_alk(gapdata = gapdata)
  
  sizecomp_stratum <- 
    calc_sizecomp_stratum(
      gapdata = gapdata, 
      cpue = cpue,
      abundance_stratum = biomass_stratum
    )
  
  sizecomp_subarea <- 
    calc_sizecomp_subarea(
      gapdata = gapdata,
      sizecomp_stratum = sizecomp_stratum
    )
  
  agecomp_stratum <-
    calc_agecomp_stratum(
      gapdata = gapdata,
      alk = alk,
      sizecomp_stratum = sizecomp_stratum
    )
  
  agecomp_region <-
    calc_agecomp_region(
      gapdata = gapdata,
      agecomp_stratum = agecomp_stratum
    )
  
  return(
    list(
      cpue = cpue, 
      biomass_stratum = biomass_stratum,
      biomass_subarea = biomass_subarea,
      alk = alk,
      sizecomp_stratum = sizecomp_stratum,
      sizecomp_subarea = sizecomp_subarea,
      agecomp_stratum = agecomp_stratum,
      agecomp_region = agecomp_region
    )
  )
  
}