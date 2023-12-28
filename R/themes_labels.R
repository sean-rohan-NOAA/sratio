#' Length label for a species
#' 
#' @param x RACEBASE species code as a numeric vector (e.g. 21740 for pollock)
#' @noRd

species_code_label <- function(x, type = "axis_label") {
  
  species_info <- data.frame(SPECIES_CODE = c(21740, 21720, 10210, 10261, 10110, 10112, 10115, 10130, 10285, 471, 68560, 68580, 69322),
                         LABEL = c(rep("Fork length (cm)", 9), "Total length (cm)", rep("Carapace width (mm)", 2), "Carapace length (mm)"),
                         COMMON_NAME = c("walleye pollock", "Pacific cod", "yellowfin sole", "northern rock sole", "arrowooth flounder", "Kamchatka flounder", "Greenland turbot", "flathead sole", "Alaska plaice", "Alaska skate", "Tanner crab", "snow crab", "red king crab"),
                         SIZE_BIN_WIDTH = c(4, 4, 3, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4))
  
  if(type == "axis_label") {
    return(species_info$LABEL[match(x, species_info$SPECIES_CODE)])
  }
   
  if(type == "common_name") {
    return(species_info$COMMON_NAME[match(x, species_info$SPECIES_CODE)])
  }
  
  if(type == "bin_width") {
    return(species_info$SIZE_BIN_WIDTH[match(x, species_info$SPECIES_CODE)])
  }
  
}


#' Length label for a species
#' 
#' @param x Numeric vector of lengths
#' @param species_code RACEBASE species code as a numeric vector (e.g. 21740 for pollock)
#' @noRd

make_size_bins <- function(x, species_code) {
  
  bin_width <- species_code_label(species_code, type = "bin_width")
  
  len_min <- min(x)
  len_max <- max(x)
  len_min <- len_min - len_min %% bin_width - bin_width
  len_max <- len_max + (bin_width - len_max %% bin_width) + bin_width
  len_breaks <- seq(len_min, len_max, bin_width)
  len_mid <- (len_breaks + bin_width/2)[1:(length(len_breaks) - 1)]
  
  out <- as.numeric(as.character(cut(x, 
                                     len_breaks, 
                                     labels = len_mid, 
                                     right = FALSE)))
  
  return(out)
  
}

