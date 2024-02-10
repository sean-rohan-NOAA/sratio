# Selectivity and catch rate analysis for 30 minute versus 15 minute hauls in the eastern Bering Sea
# Selectivity ratio, catch-at-length, and total catch rate models.
library(sratio)

# Vector of cruises to include in analyses  ----
use_cruises <- c(199501, 199801, 202101, 202201, 202301)
# use_cruises <- c(202101, 202201, 202301)
# use_cruises <- c(199801, 202101, 202201, 202301)

start_time <- Sys.time()
# 1. Retrieve groundfish and crab data from 15/30 hauls ----
# Assigns matchups
# Data will be saved as .rda files. Build and install the package to include the data as a built-in data set.
# Input: None
# Output: 
#     (1-1) ./data/data_1530.rda (built-in data set; sratio::data_1530)
#     (1-2) ./plots/sample_sizes_1530.csv (sample size table)
#     (1-3) ./plots/n_hauls.csv (hauls by year table)
#     (1-4) ./plots/sample_sizes_no_filter_1530.csv (sample size table without sample size restrictions in a haul)
source(here::here("analysis", "01_get_data.R"))


# 2. Format the built-in data ----
# Setup data for selectivity ratio and catch-at-size models 
# Input: 
#     (1-1) sratio::data_1530
# Output:
#     (2-1) ./output/catch_at_length_1530.rds (formatted for catch-at-size models)
#     (2-2) ./output/n_by_treatment_1530.rds (formatted for selectivity ratio models)
source(here::here("analysis", "02_prepare_data.R"))
  

# 3. Two-stage bootstrap samples ----
# Draw two-stage bootstrap samples for selectivity ratio and catch-at-length models using sratio::two_stage_bootstrap and sratio::nested_bootstrap
# In the two stage bootstrap, matchups are randomly drawn first then lengths are randomly drawn from each match-up and treatment. 
# Input: 
#     (2-1) ./output/catch_at_length_1530.rds
# Outputs: 
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
source(here::here("analysis", "03_bootstrap_samples.R"))


# 4. Select best binomial and beta selectivity ratio models (parallel processing) ----
# Fit GAMMs to catch comparison rate data and select the best model based on leave-one-out cross validation.
# Input: 
#     (2-2) ./output/n_by_treatment_1530.rds
# Output: 
#     (4-1) ./output/sratio_model_rmse.csv (RMSE table showing the best model for each species)
source(here::here("analysis", "04_best_selectivity_ratio_model.R"))


# 5. Bootstrap selectivity ratio models (parallel processing) ----
# For each species, fit the best model to bootstrapped sample data to estimate confidence intervals.
# Inputs:
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
#     (4-1) ./output/sratio_model_rmse.csv(RMSE table)
# Outputs:
#     (5-1) ./output/{species_code}/sratio_bootstrap_results_{species_code}.rds (bootstrap fits)
source(here::here("analysis", "05_bootstrap_selectivity_ratio_model.R"))


# 6. Plot selectivity ratio bootstrap results ----
# Plot bootstrap fits
# Inputs:
#     (5-1) ./output/{species_code}/sratio_bootstrap_results_{species_code}.rds (bootstrap fits)
# Outputs:
#     (6-1) ./plots/{species_code}_trawl_height_two_panel_ratios_n.png (catch comparison rate/selectivity ratio plot)
source(here::here("analysis", "06_plot_selectivity_ratio_model.R"))



# 7. Select best catch-at-size models (parallel processing) ----
# Fit Poisson, Negative Binomial, and Tweedie GAMMs to catch-at-size (in numbers) data and select the best models based on leave-one-out-cross validation.
# Inputs:
#     (2-1) ./output/catch_at_length_1530.rds (formatted for catch-at-length models)
# Output:
#     (7-1) ./output/cal_model_rmse.csv (RMSE table showing the best model for each species)
source(here::here("analysis", "07_best_catch_at_length_model.R"))


# 8. Bootstrap catch-at-size models (parallel processing) ----
# For each species, fit the best models to bootstrapped sample data to estimate confidence intervals.
# Inputs:
#     (3-1) ./output/{species_code}/bootstrap_samples_{species_code}.rds
#     (7-1) ./output/cal_model_rmse.csv (RMSE table showing the best model for each species)
# Outputs:
#     (8-1) ./output/{species_code}/cal_model_bootstrap_results_{species_code}.rds
source(here::here("analysis", "08_bootstrap_catch_at_length_model.R"))


# 9. Plot catch-at-size bootstrap results ----
# Plot catch-at-size fits
# Inputs: 
#     (8-1) ./output/{species_code}/cal_model_bootstrap_results_{species_code}.rds
# Outputs:
#     (9-1) ./plots/{species_code}_total_catch_gam_ratio_ribbon.png
#     (9-2) ./plots/{species_code}_total_catch_gam_ratio_lines.png
source(here::here("analysis", "09_plot_catch_at_length_model.R"))

stop_time <- Sys.time()

stop_time-start_time

# 10. Total catch model ----
# Bayesian zero-intercept linear regression model between log10(CPUE) from 15 and 30 minute hauls
# to evaluate whether there is a one-to-one relationship between CPUE from 15 and 30 minute hauls (i.e. slope = 1)
# Inputs:
#     (1-1) ./data/data_1530.rda (built-in data set; sratio::data_1530)
# Outputs:
#     (10-1) ./plots/cpue_model_density_plot.png (density plots of regression slope 95% credible intervals)
#     (10-2) ./plots/cpue_model_violin_plot.png (violin plots of regression slope  95% credible intervals)
#     (10-3) ./plots/cpue_model_boxplot.png (boxplot of regression slope 95% credible intervals)
#     (10-4) ./plots/cpue_log_model_scatterplot.png (regression fits between log10(CPUE15)~log10(CPUE30))
source(here::here("analysis/10_cpue_model.R"))

# 11. Mean bias and other metrics ---- 
#  Compare CPUE between 15 and 30 minute tows
# Inputs:
#     (1-1) ./data/data_1530.rda (built-in data set; sratio::data_1530)
# Outputs:
#     (11-1) ./plots/bias_table.csv (table of bias, RMSE, MAE by species)
source(here::here("analysis", "11_performance_metrics.R"))


# 90. Map of annual samples ----
# Inputs:
#     (1-1) ./data/data_1530.rda (built-in data set; sratio::data_1530)
# Outputs:
#     (90-1) ./plots/map_samples_by_stratum.png (samples by stratum for project plan and presentations)
source(here::here("analysis", "90_sample_map.R"))
