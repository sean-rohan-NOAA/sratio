# Effect of losing stations on crab products
# 
# This analysis examines the impact of excluding a proportion of survey stations in the EBS
# on design-based indices of abundance. Stations are excluded proportionally from each groundfish 
# stratum.
#
# Created by Sean Rohan
# March 25, 2025

library(crabpack) # 1.0.0
library(sratio)

prop_drop <- 0.25
sub_dir <- gsub(".*\\.", "", prop_drop)

crab_species_surveys <- 
  expand.grid(
    species = c("RKC", "BKC", "TANNER", "SNOW", "HYBRID", "HAIR"),
    region = c("EBS", "NBS")
  )

crab_observed_biomass <- vector(mode = "list", length = nrow(crab_species_surveys))
crab_simulated_biomass <- vector(mode = "list", length = nrow(crab_species_surveys))

for(ii in 1:nrow(crab_species_surveys)) {
  
  cat(paste0(ii, " ", crab_species_surveys$species[ii], " ", crab_species_surveys$region[ii]," ", Sys.time(), "\n"))
  
  # Load a list where each object in the list contains stations to drop
  station_draws <- readRDS(
    file = here::here("analysis", "effort_reduction", "output", 
                      sub_dir, paste0(crab_species_surveys$region[ii], "_removed_stations.rds"))
  )
  
  # Set years
  survey_years <- if(crab_species_surveys$region[ii] == "EBS") {1987:2024} else {2010:2024}
  
  # Retrieve crab data from 1987-2024
  crab_data <- crabpack::get_specimen_data(
    species = crab_species_surveys$species[ii], 
    region = crab_species_surveys$region[ii],
    years = survey_years,
    channel = "API"
  )
  
  # Calculate observed biomass and abundance
  crab_observed_biomass[[ii]] <- 
    crabpack::calc_bioabund(
      crab_data = crab_data, 
      species = crab_species_surveys$species[ii],
      region = crab_species_surveys$region[ii],
      crab_category = "all_categories"
    )
  
  simulated_reduced <- vector(mode = "list", length = length(station_draws))
  
  for(jj in 1:length(station_draws)) {
    # Drop stations ----
    sel_crab_data <- crab_data
    
    sel_crab_data$haul <- sel_crab_data$haul[!(sel_crab_data$haul$STATION_ID %in% station_draws[[jj]]), ]
    
    sel_crab_data$specimen <- sel_crab_data$specimen[!(sel_crab_data$specimen$STATION_ID %in% station_draws[[jj]]), ]
    
    # Calculate biomass and abundance after dropping stations
    simulated_reduced[[jj]] <- 
      crabpack::calc_bioabund(
        crab_data = sel_crab_data, 
        species = crab_species_surveys$species[ii],
        region = crab_species_surveys$region[ii],
        crab_category = "all_categories"
      ) |>
      dplyr::mutate(iter = jj)
  }
  
  crab_simulated_biomass[[ii]] <- do.call(rbind, simulated_reduced)
  
}

crab_biomass_results <- do.call(rbind, crab_simulated_biomass)
crab_observed_results <- do.call(rbind, crab_observed_biomass)

saveRDS(
  object = crab_biomass_results, 
  file = here::here("analysis", 
                    "effort_reduction", 
                    "output", 
                    sub_dir, 
                    paste0(
                      paste(
                        unique(crab_species_surveys$region), 
                        collapse = "_"), 
                      "_crab_biomass_district.rds")
  ),
  compress = "xz"
)

saveRDS(
  object = crab_observed_results, 
  file = here::here("analysis", 
                    "effort_reduction", 
                    "output", 
                    sub_dir, 
                    paste0(
                      paste(
                        unique(crab_species_surveys$region), 
                        collapse = "_"
                      ), "_crab_biomass_observed.rds")
  ),
  compress = "xz"
)


# Crab plots ---

crab_biomass_results <- readRDS(
  here::here("analysis", 
             "effort_reduction", 
             "output", 
             "25",
             "EBS_NBS_crab_biomass_district.rds")
)


crab_observed_results <- readRDS(
  here::here("analysis", 
             "effort_reduction", 
             "output", 
             "25",
             "EBS_NBS_crab_biomass_observed.rds")
)


# Function to make plots for each species, region, and district...
make_crab_plots <- function(reduced, observed, prop_drop, sub_dir) {


  
observed <- dplyr::select(
  observed, 
  SPECIES, 
  YEAR, 
  REGION, 
  DISTRICT, 
  CATEGORY, 
  FULL_ABUNDANCE = ABUNDANCE,
  FULL_ABUNDANCE_CV = ABUNDANCE_CV, 
  FULL_ABUNDANCE_CI = ABUNDANCE_CI, 
  FULL_BIOMASS_MT = BIOMASS_MT, 
  FULL_BIOMASS_MT_CV = BIOMASS_MT_CV,
  FULL_BIOMASS_MT_CI = BIOMASS_MT_CI, 
  FULL_BIOMASS_LBS = BIOMASS_LBS, 
  FULL_BIOMASS_LBS_CV = BIOMASS_LBS_CV
)

reduced_vs_observed <- 
  dplyr::inner_join(
    reduced, 
    observed, 
    by = c(
      "SPECIES", 
      "YEAR", 
      "REGION", 
      "DISTRICT", 
      "CATEGORY"
    )
  ) |>
  dplyr::mutate(
    PCT_ABUNDANCE = (FULL_ABUNDANCE - ABUNDANCE)/FULL_ABUNDANCE*100,
    ABS_PCT_ABUNDANCE = abs((FULL_ABUNDANCE - ABUNDANCE)/FULL_ABUNDANCE)*100,
    PCT_ABUNDANCE_CV = (FULL_ABUNDANCE_CV - ABUNDANCE_CV)/FULL_ABUNDANCE_CV*100,
    DIFF_ABUNDANCE_CV = ABUNDANCE_CV - FULL_ABUNDANCE_CV,
    PCT_BIOMASS_MT = (FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_MT = abs((FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT)*100,
    PCT_BIOMASS_MT_CV = (FULL_BIOMASS_MT_CV - BIOMASS_MT_CV)/FULL_BIOMASS_MT_CV*100,
    DIFF_BIOMASS_MT_CV = BIOMASS_MT_CV - FULL_BIOMASS_MT_CV,
    PCT_BIOMASS_LBS = (FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_LBS = abs((FULL_BIOMASS_MT - BIOMASS_MT)/FULL_BIOMASS_MT)*100,
    PCT_BIOMASS_LBS_CV = (FULL_BIOMASS_LBS_CV - BIOMASS_LBS_CV)/FULL_BIOMASS_LBS_CV*100,
    DIFF_BIOMASS_LBS_CV = BIOMASS_LBS_CV - FULL_BIOMASS_LBS_CV
  )


reduced_vs_observed_summary <- 
  reduced_vs_observed |>
  dplyr::group_by(
    SPECIES, 
    REGION, 
    DISTRICT, 
    TOTAL_AREA, 
    CATEGORY
    ) |>
  dplyr::summarise(
    MEAN_PCT_ABUNDANCE = mean(PCT_ABUNDANCE, na.rm = TRUE),
    MEAN_ABS_PCT_ABUNDANCE = mean(ABS_PCT_ABUNDANCE, na.rm = TRUE),
    MEAN_PCT_ABUNDANCE_CV = mean(PCT_ABUNDANCE_CV, na.rm = TRUE),
    MEAN_DIFF_ABUNDANCE_CV = mean(DIFF_ABUNDANCE_CV, na.rm = TRUE),
    MEAN_PCT_BIOMASS_MT = mean(PCT_BIOMASS_MT, na.rm = TRUE),
    MEAN_ABS_PCT_BIOMASS_MT = mean(ABS_PCT_BIOMASS_MT, na.rm = TRUE),
    MEAN_PCT_BIOMASS_MT_CV = mean(PCT_BIOMASS_MT_CV, na.rm = TRUE),
    MEAN_DIFF_BIOMASS_MT_CV = mean(DIFF_BIOMASS_MT_CV, na.rm = TRUE),
    MEAN_PCT_BIOMASS_LBS = mean(PCT_BIOMASS_LBS, na.rm = TRUE),
    MEAN_ABS_PCT_BIOMASS_LBS = mean(ABS_PCT_BIOMASS_LBS, na.rm = TRUE),
    MEAN_PCT_BIOMASS_LBS_CV = mean(PCT_BIOMASS_LBS_CV, na.rm = TRUE),
    MEAN_DIFF_BIOMASS_LBS_CV = mean(DIFF_BIOMASS_LBS_CV, na.rm = TRUE),
    .groups = "keep"
  )

# Plot abundance time series ----
p_abundance <- 
  ggplot() +
  geom_line(data = reduced,
            mapping = aes(
              x = YEAR, 
              y = ABUNDANCE, 
              color = paste0("Drop ", prop_drop*100, "%"), 
              group = factor(iter)),
            alpha = 0.3) +
  geom_line(data = observed,
            mapping =
              aes(
                x = YEAR,
                y = FULL_ABUNDANCE,
                color = "Observed"),
            linewidth = 1.01
  ) +
  scale_color_manual(values = c("grey50", "red")) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Abundance") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY, scales = "free")

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir,
                   "_abundance_timeseries.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_abundance)
dev.off()

# Plot abundance percent change ----
p_pct_abundance <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = PCT_ABUNDANCE),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = dplyr::group_by(reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY) |>
              dplyr::summarise(MEAN_PCT_ABUNDANCE = mean(PCT_ABUNDANCE, na.rm = TRUE)),
            mapping = aes(
              x = YEAR,
              y = MEAN_PCT_ABUNDANCE),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_PCT_ABUNDANCE),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Abundance difference (%)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_abundance_pct.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_pct_abundance)
dev.off()

# Plot absolute abundance percent change ----
p_abs_pct_abundance <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = ABS_PCT_ABUNDANCE),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = dplyr::group_by(reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY) |>
              dplyr::summarise(MEAN_ABS_PCT_ABUNDANCE = mean(ABS_PCT_ABUNDANCE, na.rm = TRUE)),
            mapping = aes(
              x = YEAR,
              y = MEAN_ABS_PCT_ABUNDANCE),
  ) +
  geom_hline(data = reduced_vs_observed_summary,
             mapping = aes(yintercept = MEAN_ABS_PCT_ABUNDANCE),
             linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Absolute abundance difference (%)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_abundance_abs_pct.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_abs_pct_abundance)
dev.off()

# Plot abundance CV time series ----
p_abundance_cv <- 
  ggplot() +
  geom_line(data = reduced,
            mapping = aes(
              x = YEAR, 
              y = ABUNDANCE_CV, 
              color = paste0("Drop ", prop_drop*100, "%"), 
              group = factor(iter)),
            alpha = 0.3) +
  geom_line(data = observed,
            mapping =
              aes(
                x = YEAR,
                y = FULL_ABUNDANCE_CV,
                color = "Observed"),
            linewidth = 1.01
  ) +
  scale_color_manual(values = c("grey50", "red")) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Abundance CV") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY, scales = "free")

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_abundance_cv_timeseries.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_abundance_cv)
dev.off()

# Plot abundance CV difference ----
p_diff_abundance_cv <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = DIFF_ABUNDANCE_CV),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = dplyr::group_by(reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY) |>
              dplyr::summarise(MEAN_DIFF_ABUNDANCE_CV = mean(DIFF_ABUNDANCE_CV, na.rm = TRUE)),
            mapping = aes(
              x = YEAR,
              y = MEAN_DIFF_ABUNDANCE_CV),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_DIFF_ABUNDANCE_CV),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = expression('Mean '*CV[reduced]-CV[full])) +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_abundance_cv_diff.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_diff_abundance_cv)
dev.off()


# Plot abundance CV percent difference ----
p_pct_abundance_cv <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = PCT_ABUNDANCE_CV),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = dplyr::group_by(reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY) |>
              dplyr::summarise(MEAN_PCT_ABUNDANCE_CV = mean(PCT_ABUNDANCE_CV, na.rm = TRUE)),
            mapping = aes(
              x = YEAR,
              y = MEAN_PCT_ABUNDANCE_CV),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_PCT_ABUNDANCE_CV),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Abundance CV difference (%)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_abundance_cv_pct.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_pct_abundance_cv)
dev.off()

# Plot biomass time series ----
p_biomass_mt <- 
  ggplot() +
  geom_line(data = reduced,
            mapping = aes(
              x = YEAR, 
              y = BIOMASS_MT, 
              color = paste0("Drop ", prop_drop*100, "%"), 
              group = factor(iter)),
            alpha = 0.3) +
  geom_line(data = observed,
            mapping =
              aes(
                x = YEAR,
                y = FULL_BIOMASS_MT,
                color = "Observed"),
            linewidth = 1.01
  ) +
  scale_color_manual(values = c("grey50", "red")) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Biomass (mt)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY, scales = "free")

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir,
                   "_biomass_mt_timeseries.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_biomass_mt)
dev.off()

# Plot biomass percent change ----
p_pct_biomass_mt <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = PCT_BIOMASS_MT),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = 
              dplyr::group_by(
                reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY
              ) |>
              dplyr::summarise(
                MEAN_PCT_BIOMASS_MT = mean(PCT_BIOMASS_MT, na.rm = TRUE)
              ),
            mapping = aes(
              x = YEAR,
              y = MEAN_PCT_BIOMASS_MT),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_PCT_BIOMASS_MT),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Biomass difference (%)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_biomass_mt_pct.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_pct_biomass_mt)
dev.off()

# Plot absolute biomass percent change ----
p_abs_pct_biomass_mt <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = ABS_PCT_BIOMASS_MT),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = 
              dplyr::group_by(
                reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY
              ) |>
              dplyr::summarise(
                MEAN_ABS_PCT_BIOMASS_MT = mean(ABS_PCT_BIOMASS_MT, na.rm = TRUE)
              ),
            mapping = aes(
              x = YEAR,
              y = MEAN_ABS_PCT_BIOMASS_MT),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_ABS_PCT_BIOMASS_MT),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Absolute biomass difference (%)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_biomass_mt_abs_pct.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_abs_pct_biomass_mt)
dev.off()

# Plot biomass CV time series ----
p_biomass_mt_cv <- 
  ggplot() +
  geom_line(data = reduced,
            mapping = aes(
              x = YEAR, 
              y = BIOMASS_MT_CV, 
              color = paste0("Drop ", prop_drop*100, "%"), 
              group = factor(iter)),
            alpha = 0.3) +
  geom_line(data = observed,
            mapping =
              aes(
                x = YEAR,
                y = FULL_BIOMASS_MT_CV,
                color = "Observed"),
            linewidth = 1.01
  ) +
  scale_color_manual(values = c("grey50", "red")) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Biomass CV") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY, scales = "free")

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_biomass_mt_cv_timeseries.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_biomass_mt_cv)
dev.off()

# Plot biomass CV difference ----
p_diff_biomass_mt_cv <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = DIFF_BIOMASS_MT_CV),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = dplyr::group_by(reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY) |>
              dplyr::summarise(MEAN_DIFF_BIOMASS_MT_CV = mean(DIFF_BIOMASS_MT_CV, na.rm = TRUE)),
            mapping = aes(
              x = YEAR,
              y = MEAN_DIFF_BIOMASS_MT_CV),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_DIFF_BIOMASS_MT_CV),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = expression('Mean '*CV[reduced]-CV[full])) +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_biomass_mt_cv_diff.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_diff_biomass_mt_cv)
dev.off()


# Plot biomass CV percent difference ----
p_pct_biomass_mt_cv <- 
  ggplot() +
  geom_jitter(
    data = reduced_vs_observed,
    mapping = aes(x = YEAR, y = PCT_BIOMASS_MT_CV),
    size = rel(0.1),
    color = "grey50",
    alpha = 0.2, 
    width = 0.25,
    height = 0) +
  geom_line(data = dplyr::group_by(reduced_vs_observed, YEAR, SPECIES, DISTRICT, CATEGORY) |>
              dplyr::summarise(MEAN_PCT_BIOMASS_MT_CV = mean(PCT_BIOMASS_MT_CV, na.rm = TRUE)),
            mapping = aes(
              x = YEAR,
              y = MEAN_PCT_BIOMASS_MT_CV),
  ) +
  # geom_hline(data = reduced_vs_observed_summary,
  #            mapping = aes(yintercept = MEAN_PCT_BIOMASS_MT_CV),
  #            linetype = 2) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Biomass CV difference (%)") +
  ggtitle(label = paste0(reduced$REGION[1], " ", reduced$SPECIES[1])) +
  theme_light() +
  theme(legend.title = element_blank()) +
  facet_grid(DISTRICT~CATEGORY)

png(filename = 
      here::here("analysis", 
                 "effort_reduction", 
                 "plots", 
                 sub_dir, 
                 reduced$REGION[1],
                 paste0(
                   reduced$REGION[1],
                   "_",
                   reduced$SPECIES[1], 
                   "_", 
                   sub_dir, 
                   "_biomass_mt_cv_pct.png"
                 )
      ),
    width = 8, 
    height = 6, 
    units = "in", 
    res = 300)
print(p_pct_biomass_mt_cv)
dev.off()

# Outputs ----
return(
  list(
  p_abundance = p_abundance,
  p_pct_abundance = p_pct_abundance,
  p_abs_pct_abundance = p_abs_pct_abundance,
  p_abundance_cv = p_abundance_cv,
  p_diff_abundance_cv = p_diff_abundance_cv,
  p_pct_abundance_cv = p_pct_abundance_cv,
  p_biomass_mt = p_biomass_mt,
  p_pct_biomass_mt = p_pct_biomass_mt,
  p_abs_pct_biomass_mt = p_abs_pct_biomass_mt,
  p_biomass_mt_cv = p_biomass_mt_cv,
  p_diff_biomass_mt_cv = p_diff_biomass_mt_cv,
  p_pct_biomass_mt_cv = p_pct_biomass_mt_cv


))

}

for(jj in 1:nrow(crab_species_surveys)) {
  
  out <- make_crab_plots(
    reduced = crab_biomass_results |>
      dplyr::filter(SPECIES == crab_species_surveys$species[[jj]],
                    REGION == crab_species_surveys$region[[jj]]), 
    observed = crab_observed_results |>
      dplyr::filter(SPECIES == crab_species_surveys$species[[jj]],
                    REGION == crab_species_surveys$region[[jj]]),
    prop_drop = prop_drop,
    sub_dir = sub_dir
  )
  
}
