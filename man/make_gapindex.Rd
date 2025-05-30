% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_gapindex.R
\name{make_gapindex}
\alias{make_gapindex}
\title{Generate gapindex results}
\usage{
make_gapindex(gapdata)
}
\arguments{
\item{gapdata}{A data frame containing survey data required for index calculations.}
}
\value{
A list containing the following data tables:
  - `cpue`: Catch per unit effort (CPUE) estimates.
  - `biomass_stratum`: Biomass estimates at the stratum level.
  - `biomass_subarea`: Biomass estimates at the subarea level.
  - `alk`: Age-length key (ALK).
  - `sizecomp_stratum`: Size composition estimates at the stratum level.
  - `sizecomp_subarea`: Size composition estimates at the subarea level.
  - `agecomp_stratum`: Age composition estimates at the stratum level.
  - `agecomp_region`: Age composition estimates at the regional level.
}
\description{
This function is a wrapper function that uses the output of `gapindex::get_data()` to generate selected data product outputs (see below).
}
\examples{
\dontrun{
library(gapindex)

gapdata <- 
gapindex::get_data(
     year_set = 1987:2024, 
     survey_set = "EBS", 
     spp_codes = 21720, 
     pull_lengths = TRUE
     )

result <- make_gapindex(gapdata)
}
}
