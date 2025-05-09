#' Length labels and size breaks for each species
#' 
#' @param x RACEBASE species code as a numeric vector (e.g. 21740 for pollock)
#' @param type Type of value to return, either "axis_label", "common_name" or "bin_width"
#' @param make_factor For common_name, should the common name be returned as an ordered factor
#' @export

species_code_label <- function(x, type = "axis_label", make_factor = FALSE) {
  
  species_info <- data.frame(SPECIES_CODE = c(21740, 21720, 10110, 10112, 10115, 10130, 10210, 10261, 10285, 471, 69322, 693221, 693222, 68580, 685801, 685802, 68560, 685601, 685602),
                         LABEL = c(rep("Fork length (cm)", 9), "Total length (cm)", rep("Carapace length (mm)", 3), rep("Carapace width (mm)", 6)),
                         COMMON_NAME = c("walleye pollock", 
                                         "Pacific cod", 
                                         "arrowtooth flounder", 
                                         "Kamchatka flounder", 
                                         "Greenland turbot", 
                                         "flathead sole", 
                                         "yellowfin sole", 
                                         "northern rock sole",  
                                         "Alaska plaice", 
                                         "Alaska skate", 
                                         "red king crab", 
                                         "red king crab (male)", 
                                         "red king crab (female)", 
                                         "snow crab",
                                         "snow crab (male)",
                                         "snow crab (female)",
                                         "Tanner crab", 
                                         "Tanner crab (male)",
                                         "Tanner crab (female)"),
                         # SIZE_BIN_WIDTH = c(6, 6, 6, 6, 6, 6, 6, 6, 6, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20),
                         SIZE_BIN_WIDTH = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8),
                         MIN_SAMPLE_SIZE = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10)
                         )
  
  if(x[1] == "all") {
    return(species_info)
  }
  
  if(type == "axis_label") {
    return(species_info$LABEL[match(x, species_info$SPECIES_CODE)])
  }
   
  if(type == "common_name") {
    
    out <- species_info$COMMON_NAME[match(x, species_info$SPECIES_CODE)]
    
    if(make_factor) {
      
      out <- factor(out, levels = species_info$COMMON_NAME)
      
    }
    
    return(out)
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