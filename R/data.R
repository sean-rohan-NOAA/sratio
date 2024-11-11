#' 15/30 Comparison Data
#' 
#' List containing haul, catch, length, and crab size data for the 15/30 comparison project.
#' 
#' @docType data
#' @format An object of class 'list'
#' 
"data_1530"


#' gapindex data for the 15/30 comparison
#' 
#' List containing RACEBASE tables created using the gapindex package for the 15/30 comparison project.
#' 
#' @docType data
#' @format An object of class 'list'
#' 
"gapindex_1530"


#' Shelf/slope Comparison Data
#' 
#' List containing haul, catch, and size (fish length, crab carapace width, crab carapace length) for the shelf/slope comparison project.
#' 
#' @docType data
#' @format An object of class 'list'
#' 
"data_ss"


#' Crab size from 15/30 minute tows in 1995 and 1998
#' 
#' Crab size measurement data from paired 15 and 30 minute tows conducted in the eastern Bering Sea in 1995 and 1998 aboard the FV Arcturus and FV Aldebaran. Sampling protocols and methods are described in Goddard et al. (1997) and Somerton et al. (2002).
#' 
#' @docType data
#' @usage data(crab_size_1995_1998)
#' @format A data frame with 12 variables:
#'   \describe{
#'     \item{HAULJOIN}{A unique identifier for each CRUISE/VESSEL/HAUL combination.}
#'     \item{CRUISE}{Six-digit number identifying the cruise, where the first four digits are the year of the cruise and last four digits are region-specific sequential codes.}
#'     \item{VESSEL}{Unique numeric vessel code.}
#'     \item{HAUL}{Sequential haul number identifying a unique haul for a VESSEL/CRUISE combination}
#'     \item{SPECIES_CODE}{Numeric species code}
#'     \item{SEX}{Sex, where 1 = male, 2 = female, 3 = unknown/unspecified}
#'     \item{LENGTH}{Carapace length in millimeters}
#'     \item{WIDTH}{Carapace width in millimeters}
#'     \item{SHELL_CONDITION}{Shell condition}
#'     \item{SAMPLING_FACTOR}{Sampling factor}
#'     \item{CHELA_HEIGHT}{Chela height in millimeters}
#'     \item{STATIONID}{Survey station name}
#'   }
#' @references Goddard, P.D., 1997. The effects of tow duration and subsampling on CPUE, species composition and length distributions of bottom trawl survey catches. Master of Science Thesis, University of Washington, Seattle, WA, 119 pp.
#'  Somerton, D.A., Otto, R.S., and Syrjala, S.E. 2002. Can changes in tow duration on bottom trawl surveys lead to changes in CPUE and mean size? Fisheries Research, 55(1–3), 63–70. https://doi.org/10.1016/S0165-7836(01)00293-4
#' 
'crab_size_1995_1998'


#' Blocks of 15/30 minute tows from 1998
#' 
#' Key identifying blocks of 15/30 minute tows adapted from Bob Otto's tow pair spreadsheet (2000).
#' 
#' @docType data
#' @usage data(otto_key_1998)
#' @format A data frame with 6 variables:
#'   \describe{
#'     \item{TOW_PAIR}{Value identifying the block to which a tow was assigned.}
#'     \item{EXPERIMENTAL_TOW_TYPE}{Numerical tow type code.}
#'     \item{CRUISE}{Six-digit number identifying the cruise, where the first four digits are the year of the cruise and last four digits are region-specific sequential codes.}
#'     \item{VESSEL}{Unique numeric vessel code.}
#'     \item{HAUL}{Sequential haul number identifying a unique haul for a VESSEL/CRUISE combination}
#'     \item{TOW_TYPE}{Description of numberical tow type}
#'   }
'otto_key_1998'

#' Pacific cod selectivity ratios
#' 
#' Pacific cod selectivity ratio values calculated from 15/30 tow duration hauls.
#' 
#' @docType data
#' @format An object of class 'data.frame'
#'   \describe{
#'     \item{species_code}{Numeric species code for Pacific cod (21720)}
#'     \item{block}{Tow pair}
#'     \item{size}{Size bin; fork length in cm}
#'     \item{count1}{Number of fish lengthed in the 30 minute haul.}
#'     \item{count2}{Number of fish lengthed in the 15 minute haul.}
#'     \item{total_count}{Total number of fish from both hauls conbined; count1 + count2}
#'     \item{effort1}{Area swept in square kilometers for the 30 minute haul.}
#'     \item{effort2}{Area swept in square kilometes for the 15 minute haul.}
#'     \item{sampling_factor1}{Multiplier to scale count to total catch for the 30 minute haul, i.e., total number of fish at size = count * sampling_factor. Inverse of the estimated sampled propoprtion.}
#'     \item{sampling_factor2}{Multiplier to scale count to total catch for the 15 minute haul, i.e., total number of fish at size = count * sampling_factor Inverse of the estimated sampled propoprtion.}
#'     \item{p}{Catch comparison rate estimated using sratio::selectivity_ratio}
#'     \item{s}{Selectivity ratio between the 30 and 15 minute haul estimated using sratio::selectivity_ratio}
#'   }
"cod_sratio"