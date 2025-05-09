library(akgfmaps)
library(sratio)
library(crabpack)
library(gapindex)
library(cowplot)

# Control pars
prop_drop <- 0.5
sub_dir <- gsub(".*\\.", "", format(prop_drop, nsmall = 2))
n_iter <- 100
seed <- 1337

# Connect 
channel <- gapindex::get_connected()

survey_opts <- c("EBS", "NBS") # EBS or NBS

for(vv in 1:length(survey_opts)) {

    survey_set <- survey_opts[vv]
  
  # Setup output directories
  dir.create(path = here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set), 
             showWarnings = FALSE, 
             recursive = TRUE)
  dir.create(path = here::here("analysis", "effort_reduction", "output", sub_dir), 
             showWarnings = FALSE, 
             recursive = TRUE)
  
  survey_years <- if(survey_set == "EBS") {1987:2024} else {2010:2024}
  survey_definition_id <- ifelse(survey_set == "EBS", 98, 143)
  

  plot_area_ids <- if(survey_set == "EBS") {c(99900, 99901)} else {99902}
  loop_area_names <- if(survey_set == "EBS") {c("EBS Standard Plus NW", "EBS Standard")} else{"NBS"}
  
  source(here::here("analysis", "effort_reduction", "201_effort_reduction_groundfish.R"))
  source(here::here("analysis", "effort_reduction", "202_plot_effort_reduction_groundfish.R"))
  source(here::here("analysis", "effort_reduction", "203_effort_reduction_crab.R"))
  source(here::here("analysis", "effort_reduction", "204_plot_effort_reduction_crab.R"))
  
}

source(here::here("analysis", "effort_reduction", "205_effort_reduction_levels.R"))
