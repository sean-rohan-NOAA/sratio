library(sratio)
library(cowplot)

# Control pars
prop_drop <- 0.25
sub_dir <- gsub(".*\\.", "", prop_drop)
survey_set <- "EBS" # EBS or NBS
survey_years <- if(survey_set == "EBS") {1987:2024} else {2010:2024}
survey_definition_id <- ifelse(survey_set == "EBS", 98, 143)

plot_area_ids <- if(survey_set == "EBS") {c(99900, 99901)} else {99902}
loop_area_names <- if(survey_set == "EBS") {c("EBS Standard Plus NW", "EBS Standard")} else{"NBS"}

# Load species list
akfin_species <- 
  read.csv(
    file = here::here("analysis", "effort_reduction", "data", "akfin_biomass_taxa_by_survey.csv")
  ) |>
  dplyr::filter(SURVEY_DEFINITION_ID == survey_definition_id, 
                PLOT)

# Calculate summary statistics -----

biomass_subarea_results <- 
  readRDS(
    here::here("analysis", "effort_reduction", "output", sub_dir, 
               paste0(survey_set, "_fish_biomass_subarea.rds")
    )
  )

observed <- 
  readRDS(
    here::here("analysis", "effort_reduction", "output", sub_dir, 
               paste0(survey_set, "_fish_observed.rds")
    )
  )

# Combine outputs into a data frame
reduced_sampling_results <- do.call(rbind, biomass_subarea_results)

# Calculate CVs
reduced_sampling_results$CV <- sqrt(reduced_sampling_results$BIOMASS_VAR)/reduced_sampling_results$BIOMASS_MT

reduced_pct_change <- 
  reduced_sampling_results |>
  dplyr::inner_join(
    observed$biomass_subarea |>
      dplyr::select(SPECIES_CODE, 
                    AREA_ID,
                    YEAR,
                    FULL_BIOMASS_MT = BIOMASS_MT, 
                    FULL_BIOMASS_VAR = BIOMASS_VAR,
                    FULL_CV = CV),
    by = c("YEAR", "SPECIES_CODE", "AREA_ID")
  ) |>
  dplyr::mutate(
    PCT_CV = (CV-FULL_CV)/FULL_CV*100,
    DIFF_CV = CV-FULL_CV,
    PCT_BIOMASS_MT = (BIOMASS_MT-FULL_BIOMASS_MT)/FULL_BIOMASS_MT*100,
    ABS_PCT_BIOMASS_MT = (BIOMASS_MT-FULL_BIOMASS_MT)/FULL_BIOMASS_MT*100
  ) |>
  dplyr::filter(AREA_ID %in% plot_area_ids)


reduced_pct_change_year <- 
  reduced_pct_change |>
  dplyr::group_by(SPECIES_CODE, AREA_ID, YEAR) |>
  dplyr::summarise(
    MEAN_PCT_CV = mean(PCT_CV, na.rm = TRUE),
    MEAN_DIFF_CV = mean(DIFF_CV, na.rm = TRUE),
    MEAN_ABS_PCT_BIOMASS_MT = mean(abs(PCT_BIOMASS_MT), na.rm = TRUE),
    MEAN_PCT_BIOMASS_MT = mean(PCT_BIOMASS_MT, na.rm = TRUE)
  )

reduced_pct_change_overall <- 
  reduced_pct_change |>
  dplyr::group_by(SPECIES_CODE, AREA_ID) |>
  dplyr::summarise(
    MEAN_PCT_CV = mean(PCT_CV, na.rm = TRUE),
    MEAN_DIFF_CV = mean(DIFF_CV, na.rm = TRUE),
    MEAN_ABS_PCT_BIOMASS_MT = mean(abs(PCT_BIOMASS_MT), na.rm = TRUE),
    MEAN_PCT_BIOMASS_MT = mean(PCT_BIOMASS_MT, na.rm = TRUE)
  )


# Loop through Standard and Standard Plus NW 
loop_area_id <- plot_area_ids

for(kk in 1:length(loop_area_id)) {
  
  # Plot: Biomass timeseries ----
  
  p_biomass <- 
    ggplot() +
    geom_line(data = dplyr::filter(
      reduced_sampling_results, 
      AREA_ID == loop_area_id[kk]) |>
        dplyr::inner_join(akfin_species),
      mapping = aes(
        x = YEAR, 
        y = BIOMASS_MT, 
        color = paste0("Drop ", prop_drop*100, "%"), 
        group = factor(iter)),
      alpha = 0.3) +
    geom_line(data = 
                dplyr::filter(
                  observed$biomass_subarea, 
                  AREA_ID == loop_area_id[kk]) |>
                dplyr::inner_join(akfin_species),
              mapping = 
                aes(
                  x = YEAR, 
                  y = BIOMASS_MT, 
                  color = "Observed"),
              linewidth = 1.02
    ) +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Biomass (mt)") +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE, scales = "free_y")
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0(survey_set, "_fish_", loop_area_id[kk], "_biomass_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_biomass)
  dev.off()
  
  # Plot: CV timeseries ----
  p_cv <- 
    ggplot() +
    geom_line(data = dplyr::filter(
      reduced_sampling_results, 
      AREA_ID == loop_area_id[kk]) |>
        dplyr::inner_join(akfin_species),
      mapping = aes(
        x = YEAR, 
        y = CV, 
        color = paste0("Drop ", prop_drop*100, "%"), 
        group = factor(iter)),
      alpha = 0.3) +
    geom_line(data = 
                dplyr::filter(
                  observed$biomass_subarea, 
                  AREA_ID == loop_area_id[kk]) |>
                dplyr::inner_join(akfin_species),
              mapping = 
                aes(
                  x = YEAR, 
                  y = CV, 
                  color = "Observed"),
              linewidth = 1.02
    ) +
    scale_color_manual(values = c("grey50", "red")) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "CV", limits = c(0, NA)) +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE, scales = "free_y")
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0(survey_set, "_fish_", loop_area_id[kk], "_cv_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_cv)
  dev.off()
  
  # Plot: Biomass percent change timeseries ----
  p_biomass_pct <-
    ggplot() +
    geom_jitter(
      data =
        dplyr::filter(
          reduced_pct_change,
          AREA_ID == loop_area_id[kk]
        ) |>
        dplyr::inner_join(akfin_species),
      mapping = aes(x = YEAR, y = PCT_BIOMASS_MT),
      size = rel(0.2),
      color = "grey50",
      alpha = 0.3, 
      width = 0.25,
      height = 0) +
    geom_line(data = 
                dplyr::filter(
                  reduced_pct_change_year,
                  AREA_ID == loop_area_id[kk]
                ) |>
                dplyr::inner_join(akfin_species),
              mapping = aes(
                x = YEAR, 
                y = MEAN_PCT_BIOMASS_MT),
              size = 1.05
    ) +
    geom_hline(data = 
                 dplyr::filter(
                   reduced_pct_change_overall,
                   AREA_ID == loop_area_id[kk]
                 ) |>
                 dplyr::inner_join(akfin_species),
               mapping = 
                 aes(yintercept = MEAN_PCT_BIOMASS_MT),
               linetype = 2) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Biomass difference (%)") +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE, scales = "free_y")
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0(survey_set, "_fish_",  loop_area_id[kk], "_biomass_pct_ts.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_biomass_pct)
  dev.off()
  
  # Plot: Biomass absolute percent change timeseries ----
  # p_biomass_abs_pct <-
  #   ggplot() +
  #   geom_jitter(
  #     data =
  #       dplyr::filter(
  #         reduced_pct_change,
  #         AREA_ID == loop_area_id[kk]
  #       ) |>
  #       dplyr::inner_join(akfin_species),
  #     mapping = aes(x = YEAR, y = ABS_PCT_BIOMASS_MT),
  #     size = rel(0.3),
  #     color = "grey50",
  #     alpha = 0.2, 
  #     width = 0.25,
  #     height = 0) +
  #   geom_line(data = 
  #               dplyr::filter(
  #                 reduced_pct_change_year,
  #                 AREA_ID == loop_area_id[kk]
  #               ) |>
  #               dplyr::inner_join(akfin_species),
  #             mapping = aes(x = YEAR, y = MEAN_ABS_PCT_BIOMASS_MT)) +
  #   geom_hline(data = 
  #                dplyr::filter(
  #                  reduced_pct_change_overall,
  #                  AREA_ID == loop_area_id[kk]
  #                ) |>
  #                dplyr::inner_join(akfin_species),
  #              mapping = 
  #                aes(yintercept = MEAN_ABS_PCT_BIOMASS_MT),
  #              linetype = 2) +
  #   scale_x_continuous(name = "Year") +
  #   scale_y_continuous(name = "Mean biomass change (%)", limits = c(0, NA)) +
  #   ggtitle(label = loop_area_names[kk]) +
  #   theme_light() +
  #   theme(legend.title = element_blank()) +
  #   facet_wrap(~SPECIES_CODE, scales = "free_y")
  # 
  # png(filename = 
  #       here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
  #                  paste0(survey_set, "_fish_",  loop_area_id[kk], "_biomass_abs_pct_ts.png")),
  #     width = 8, 
  #     height = 6, 
  #     units = "in", 
  #     res = 300)
  # print(p_biomass_abs_pct)
  # dev.off()
  
  # Plot: CV change ----
  
  p_cv_change <- 
    ggplot() +
    geom_jitter(
      data =
        dplyr::filter(
          reduced_pct_change,
          AREA_ID == loop_area_id[kk]
        ) |>
        dplyr::inner_join(akfin_species),
      mapping = aes(x = YEAR, y = DIFF_CV),
      size = rel(0.2),
      color = "grey50",
      alpha = 0.5) +
    geom_line(
      data = 
        dplyr::filter(
          reduced_pct_change_year,
          AREA_ID == loop_area_id[kk]
        ) |>
        dplyr::inner_join(akfin_species),
      mapping = 
        aes(x = YEAR, 
            y = MEAN_DIFF_CV)
    ) +
    geom_hline(
      data = 
        dplyr::filter(
          reduced_pct_change_overall,
          AREA_ID == loop_area_id[kk]
        ) |>
        dplyr::inner_join(akfin_species),
      mapping = 
        aes(yintercept = MEAN_DIFF_CV),
      linetype = 2) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = expression('Mean '*CV[reduced]-CV[full])) +
    ggtitle(label = loop_area_names[kk]) +
    theme_light() +
    theme(legend.title = element_blank()) +
    facet_wrap(~SPECIES_CODE, scales = "free_y")
  
  png(filename = 
        here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                   paste0(survey_set, "_fish_", loop_area_id[kk], "_cv_diff.png")),
      width = 8, 
      height = 6, 
      units = "in", 
      res = 300)
  print(p_cv_change)
  dev.off()
  
  for(ll in 1:nrow(akfin_species)) {
    
    # Plot: Biomass timeseries ----
    
    p_biomass <- 
      ggplot() +
      geom_line(data = dplyr::filter(
        reduced_sampling_results, 
        AREA_ID == loop_area_id[kk],
        SPECIES_CODE == akfin_species$SPECIES_CODE[ll]),
        mapping = aes(
          x = YEAR, 
          y = BIOMASS_MT, 
          color = paste0("Drop ", prop_drop*100, "%"), 
          group = factor(iter)),
        alpha = 0.3) +
      geom_line(data = 
                  dplyr::filter(
                    observed$biomass_subarea, 
                    AREA_ID == loop_area_id[kk],
                    SPECIES_CODE == akfin_species$SPECIES_CODE[ll]),
                mapping = 
                  aes(
                    x = YEAR, 
                    y = BIOMASS_MT, 
                    color = "Observed"),
                linewidth = 1.02
      ) +
      scale_color_manual(values = c("grey50", "red")) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "Biomass (mt)") +
      theme_light() +
      theme(legend.title = element_blank(),
            legend.position = "none")
    
    # Plot: CV timeseries ----
    p_cv <- 
      ggplot() +
      geom_line(data = dplyr::filter(
        reduced_sampling_results, 
        AREA_ID == loop_area_id[kk],
        SPECIES_CODE == akfin_species$SPECIES_CODE[ll]),
        mapping = aes(
          x = YEAR, 
          y = CV, 
          color = paste0("Drop ", prop_drop*100, "%"), 
          group = factor(iter)),
        alpha = 0.3) +
      geom_line(data = 
                  dplyr::filter(
                    observed$biomass_subarea, 
                    AREA_ID == loop_area_id[kk],
                    SPECIES_CODE == akfin_species$SPECIES_CODE[ll]),
                mapping = 
                  aes(
                    x = YEAR, 
                    y = CV, 
                    color = "Observed"),
                linewidth = 1.02
      ) +
      scale_color_manual(values = c("grey50", "red")) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "CV", limits = c(0, NA)) +
      theme_light() +
      theme(legend.title = element_blank(),
            legend.position = "none")

    # Plot: Biomass percent change timeseries ----
    p_biomass_pct <-
      ggplot() +
      geom_jitter(
        data =
          dplyr::filter(
            reduced_pct_change,
            AREA_ID == loop_area_id[kk],
            SPECIES_CODE == akfin_species$SPECIES_CODE[ll]
          ),
        mapping = aes(x = YEAR, y = PCT_BIOMASS_MT),
        size = rel(0.2),
        color = "grey50",
        alpha = 0.3, 
        width = 0.25,
        height = 0) +
      geom_line(data = 
                  dplyr::filter(
                    reduced_pct_change_year,
                    AREA_ID == loop_area_id[kk],
                    SPECIES_CODE == akfin_species$SPECIES_CODE[ll]
                  ),
                mapping = aes(
                  x = YEAR, 
                  y = MEAN_PCT_BIOMASS_MT),
                size = 1.05
      ) +
      geom_hline(data = 
                   dplyr::filter(
                     reduced_pct_change_overall,
                     AREA_ID == loop_area_id[kk],
                     SPECIES_CODE == akfin_species$SPECIES_CODE[ll]
                   ),
                 mapping = 
                   aes(yintercept = MEAN_PCT_BIOMASS_MT),
                 linetype = 2) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "Biomass difference (%)") +
      theme_light() +
      theme(legend.title = element_blank())
    
    # Plot: CV change ----
    
    p_cv_change <- 
      ggplot() +
      geom_jitter(
        data =
          dplyr::filter(
            reduced_pct_change,
            AREA_ID == loop_area_id[kk],
            SPECIES_CODE == akfin_species$SPECIES_CODE[ll]
          ),
        mapping = aes(x = YEAR, y = DIFF_CV),
        size = rel(0.2),
        color = "grey50",
        alpha = 0.5) +
      geom_line(
        data = 
          dplyr::filter(
            reduced_pct_change_year,
            AREA_ID == loop_area_id[kk],
            SPECIES_CODE == akfin_species$SPECIES_CODE[ll]
          ),
        mapping = 
          aes(x = YEAR, 
              y = MEAN_DIFF_CV)
      ) +
      geom_hline(
        data = 
          dplyr::filter(
            reduced_pct_change_overall,
            AREA_ID == loop_area_id[kk],
            SPECIES_CODE == akfin_species$SPECIES_CODE[ll]
          ),
        mapping = 
          aes(yintercept = MEAN_DIFF_CV),
        linetype = 2) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = expression('Mean '*CV[reduced]-CV[full])) +
      theme_light() +
      theme(legend.title = element_blank())
    
    title <- cowplot::ggdraw() + 
      cowplot::draw_label(paste0(loop_area_names[kk], " ", akfin_species$COMMON_NAME[ll]) , fontface = 'bold', size = 14)
    
    png(filename = 
          here::here("analysis", "effort_reduction", "plots", sub_dir, survey_set, 
                     paste0(survey_set, "_fish_", loop_area_id[kk], "_", akfin_species$SPECIES_CODE[ll], "_diff.png")),
        width = 8, 
        height = 6, 
        units = "in", 
        res = 300)
    print(cowplot::plot_grid(
      plot_grid(title),
      plot_grid(
      p_biomass,
      p_biomass_pct,
      p_cv,
      p_cv_change),
      nrow = 2,
      rel_heights = c(0.1,1)
    ))
    dev.off()
    

    
  }
  
}
