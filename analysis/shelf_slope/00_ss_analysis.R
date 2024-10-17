# Selectivity and catch rate analysis for 30 minute versus 15 minute hauls in the eastern Bering Sea
# Selectivity ratio, catch-at-length, and total catch rate models.
library(sratio)

# Vector of cruises to include in analyses  ----
use_cruises <- c(202301, 202401)

# 1. Retrieve groundfish and crab data from shelf/slope comparison hauls ----
# Assigns matchups
# Data will be saved as .rda files. Build and install the package to include the data as a built-in data set.
# Input: None
# Output: 
#     (1-1) ./data/data_ss.rda (built-in data set; sratio::data_ss)
#     (1-2) ./plots/sample_sizes_ss.csv (sample size table)
#     (1-3) ./plots/n_hauls.csv (hauls by year table)
#     (1-4) ./plots/sample_sizes_no_filter_ss.csv (sample size table without sample size restrictions in a haul)
source(here::here("analysis", "shelf_slope",  "01_ss_get_data.R"))


#--------------------------------------------------------------------------------------------------#
#--------------------------- REBUILD THE PACKAGE BEFORE CONTINUING TO 2 ---------------------------#
#--------------------------------------------------------------------------------------------------#

start_time <- Sys.time()

# 2. Format the built-in data ----
# Setup data for selectivity ratio and catch-at-size models 
# Input: 
#     (1-1) sratio::data_ss
# Output:
#     (2-1) ./output/catch_at_length_ss.rds (formatted for catch-at-size models)
#     (2-2) ./output/n_by_treatment_ss.rds (formatted for selectivity ratio models)
source(here::here("analysis", "shelf_slope", "02_ss_prepare_data.R"))
  

# 3. Two-stage bootstrap samples ----
# Draw two-stage bootstrap samples for selectivity ratio and catch-at-length models using sratio::two_stage_bootstrap and sratio::nested_bootstrap
# In the two stage bootstrap, matchups are randomly drawn first then lengths are randomly drawn from each match-up and treatment. 
# Input: 
#     (2-1) ./output/catch_at_length_ss.rds
# Outputs: 
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
source(here::here("analysis", "shelf_slope", "03_ss_bootstrap_samples.R"))


# 4. Select best binomial and beta selectivity ratio models (parallelized) ----
# Fit GAMMs to catch comparison rate data and select the best model based on leave-one-out cross validation.
# Input: 
#     (2-2) ./output/n_by_treatment_ss.rds
# Output: 
#     (4-1) ./output/sratio_model_rmse.csv (RMSE table showing the best model for each species)
source(here::here("analysis", "shelf_slope", "04_ss_best_selectivity_ratio_model.R"))


# 5. Bootstrap selectivity ratio models (parallelized) ----
# For each species, fit the best model to bootstrapped sample data to estimate confidence intervals.
# Inputs:
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
#     (4-1) ./output/sratio_model_rmse.csv(RMSE table)
# Outputs:
#     (5-1) ./output/{species_code}/sratio_bootstrap_results_{species_code}.rds (bootstrap fits)
source(here::here("analysis", "shelf_slope", "05_ss_bootstrap_selectivity_ratio_model.R"))


# 6. Plot selectivity ratio bootstrap results ----
# Plot bootstrap fits
# Inputs:
#     (5-1) ./output/{species_code}/sratio_bootstrap_results_{species_code}.rds (bootstrap fits)
# Outputs:
#     (6-1) ./plots/{species_code}_sratio_two_panel.png (catch comparison rate/selectivity ratio plot)
#     (6-2) ./plots/{species_code}_sratio_three_panel.png (catch comparison rate/selectivity ratio plot w/ sample histogram)
source(here::here("analysis", "shelf_slope", "06_ss_plot_selectivity_ratio_model.R"))



# 7. Select best catch-at-size models (parallelized) ----
# Fit Poisson, Negative Binomial, and Tweedie GAMMs to catch-at-size (in numbers) data and select the best models based on leave-one-out-cross validation.
# Inputs:
#     (2-1) ./output/catch_at_length_ss.rds (formatted for catch-at-length models)
# Output:
#     (7-1) ./output/cal_model_rmse.csv (RMSE table showing the best model for each species)
source(here::here("analysis", "shelf_slope", "07_ss_best_catch_at_length_model.R"))


# 8. Bootstrap catch-at-size models (parallelized) ----
# For each species, fit the best models to bootstrapped sample data to estimate confidence intervals.
# Inputs:
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
#     (7-1) ./output/cal_model_rmse.csv (RMSE table showing the best model for each species)
# Outputs:
#     (8-1) ./output/{species_code}/sccal_model_bootstrap_results_{species_code}.rds
source(here::here("analysis", "shelf_slope", "08_ss_bootstrap_catch_at_length_model.R"))


# 9. Plot selectivity conditional on catch-at-length (SCCAL) bootstrap results ----
# Plot catch-at-size fits
# Inputs: 
#     (8-1) ./output/{species_code}/sccal_model_bootstrap_results_{species_code}.rds
# Outputs:
#     (9-1) ./plots/{species_code}_sccal_ratio.png
source(here::here("analysis", "shelf_slope", "09_ss_plot_catch_at_length_model.R"))


# 10. Plot SCCAL and SR on the same plots ----
# Inputs: 
#     (5-1) ./output/{species_code}/sratio_bootstrap_results_{species_code}.rds (bootstrap fits)
#     (8-1) ./output/{species_code}/sccal_model_bootstrap_results_{species_code}.rds
# source(here::here("analysis", "shelf_slope", "10_plot_sratio_sccal.R"))
# 
# stop_time <- Sys.time()
# 
# stop_time-start_time


# 11. Calculate performance metrics ----
# Inputs: 
#     (1-1) sratio::data_ss
#     (5-1) ./output/{species_code}/sratio_bootstrap_results_{species_code}.rds (bootstrap fits)
#     (8-1) ./output/{species_code}/sccal_model_bootstrap_results_{species_code}.rds
# Outputs:
#     (11-1) ./output/model_performance_sratio_{species_code}.rds
#     (11-2) ./output/model_performance_sccal_{species_code}.rds

source(here::here("analysis", "shelf_slope", "11_ss_calculate_performance_metrics.R"))


# 12. Plot performance metrics ----




# 13. Total catch model ----
# Bayesian zero-intercept linear regression model between log10(CPUE) from shelf and slope tows
# to evaluate whether there is a one-to-one relationship between CPUE from shelf and slope tows (i.e. slope = 1)
# Inputs:
#     (1-1) ./data/data_ss.rda (built-in data set; sratio::data_ss)
# Outputs:
#     (11-1) ./plots/cpue_model_density_plot.png (density plots of regression slope 95% credible intervals)
#     (11-2) ./plots/cpue_model_violin_plot.png (violin plots of regression slope  95% credible intervals)
#     (11-3) ./plots/cpue_model_boxplot.png (boxplot of regression slope 95% credible intervals)
#     (11-4) ./plots/cpue_log_model_scatterplot.png (regression fits between log10(CPUE15)~log10(CPUE30))
source(here::here("analysis", "shelf_slope", "11_ss_cpue_model.R"))


# 12. Mean bias and other metrics ---- 
#  Compare CPUE between shelf and slope tows
# Inputs:
#     (1-1) ./data/data_ss.rda (built-in data set; sratio::data_ss)
# Outputs:
#     (12-1) ./plots/bias_table.csv (table of bias, RMSE, MAE by species)
source(here::here("analysis", "shelf_slope", "12_ss_performance_metrics.R"))


# 90. Map of annual samples ----
# Inputs:
#     (1-1) ./data/data_ss.rda (built-in data set; sratio::data_ss)
# Outputs:
#     (90-1) ./plots/map_samples_by_stratum.png (samples by stratum for project plan and presentations)
#     (90-2) ./plots/sample_map_1995_2023.png (multi-panel sample map/samples by year)
source(here::here("analysis", "shelf_slope", "90_ss_sample_map.R"))
