# Selectivity and catch rate analysis for 30 minute versus 15 minute hauls in the eastern Bering Sea
# Selectivity ratio, catch-at-length, and total catch rate models.
library(sratio)

# Vector of cruises to include in analyses  ----
use_cruises <- c(199501, 199801, 202101, 202201, 202301, 202401)

# 1. Retrieve groundfish and crab data from 15/30 hauls ----
# Assigns matchups
# Data will be saved as .rda files. Build and install the package to include the data as a built-in data set.
# Input: None
# Output: 
#     (1-1) ./data/data_1530.rda (built-in data set; sratio::data_1530)
#     (1-2) ./plots/sample_sizes_1530.csv (sample size table)
#     (1-3) ./plots/n_hauls.csv (hauls by year table)
#     (1-4) ./plots/sample_sizes_no_filter_1530.csv (sample size table without sample size restrictions in a haul)
source(here::here("analysis", "15_30",  "01_get_data.R"))


#--------------------------------------------------------------------------------------------------#
#--------------------------- REBUILD THE PACKAGE BEFORE CONTINUING TO 2 ---------------------------#
#--------------------------------------------------------------------------------------------------#

start_time <- Sys.time()

# 2. Format the built-in data ----
# Setup data for selectivity ratio and catch-at-size models 
# Input: 
#     (1-1) sratio::data_1530
# Output:
#     (2-1) ./output/catch_at_length_1530.rds (formatted for catch-at-size models)
#     (2-2) ./output/n_by_treatment_1530.rds (formatted for selectivity ratio models)
source(here::here("analysis", "15_30", "02_prepare_data.R"))
  

# 3. Two-stage bootstrap samples ----
# Draw two-stage bootstrap samples for selectivity ratio and catch-at-length models using sratio::two_stage_bootstrap and sratio::nested_bootstrap
# In the two stage bootstrap, matchups are randomly drawn first then lengths are randomly drawn from each match-up and treatment. 
# Input: 
#     (2-1) ./output/catch_at_length_1530.rds
# Outputs: 
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
source(here::here("analysis", "15_30", "03_bootstrap_samples.R"))

# 90. Map of annual samples ----
# Inputs:
#     (1-1) ./data/data_1530.rda (built-in data set; sratio::data_1530)
# Outputs:
#     (90-1) ./plots/map_samples_by_stratum.png (samples by stratum for project plan and presentations)
#     (90-2) ./plots/sample_map_1995_2023.png (multi-panel sample map/samples by year)
source(here::here("analysis", "15_30", "90_sample_map.R"))
