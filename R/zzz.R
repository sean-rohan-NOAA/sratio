#' @useDynLib sratio
.onUnload <- function (lib) {
  library.dynam.unload("sratio", lib)
}

#' Internal function to make a wide-format CPUE data frame
#' 
#' @noRd

make_cpue_wide <- function(x) {
  
  cpue_kgkm2 <- x$catch |>
    dplyr::inner_join(x$haul) |>
    dplyr::mutate(CPUE_KGKM2 = WEIGHT/AREA_SWEPT_KM2) |>
    dplyr::select(TREATMENT, MATCHUP, CPUE_KGKM2, SPECIES_CODE) |>
    tidyr::pivot_wider(values_from = "CPUE_KGKM2", 
                       names_from = "TREATMENT", 
                       names_prefix = "CPUE_KGKM2_", 
                       values_fill = 0)
  
  cpue_number_km2 <- x$catch |>
    dplyr::inner_join(x$haul) |>
    dplyr::mutate(CPUE_NUMBER_KM2 = NUMBER_FISH/AREA_SWEPT_KM2) |>
    dplyr::select(TREATMENT, MATCHUP, CPUE_NUMBER_KM2, SPECIES_CODE) |>
    tidyr::pivot_wider(values_from = "CPUE_NUMBER_KM2", 
                       names_from = "TREATMENT", 
                       names_prefix = "CPUE_NUMBER_KM2_", 
                       values_fill = 0)
  
  area_swept_km2 <- x$catch |>
    dplyr::inner_join(x$haul) |>
    dplyr::select(TREATMENT, MATCHUP, SPECIES_CODE, AREA_SWEPT_KM2) |>
    tidyr::pivot_wider(values_from = "AREA_SWEPT_KM2", 
                       names_from = "TREATMENT", 
                       names_prefix = "AREA_SWEPT_KM2_", 
                       values_fill = 0)
  
  return(dplyr::inner_join(cpue_kgkm2, cpue_number_km2))
  
}