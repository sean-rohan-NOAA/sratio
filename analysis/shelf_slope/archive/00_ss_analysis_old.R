# Selectivity ratio and catchabiltiy calibration for shelf/slope comparison

# 1. Setup ----
library(sratio)

n_cores <- 4 # Cores for parallel processing
bin_width <- c(5, 5, 4, 4, 5, 5)  # May need to change for different species
species_codes <- c(21740, 21720, 10130, 10115, 10110, 10112, 471) #,68580,  68560, 69322)
measurement_label <- c(rep("Fork length (cm)", 6), "Total length (cm)") #, "Carapace length (mm)", "Carapace width (mm)", rep("Fork length (cm)", 2), "Carapace width (mm)")
seed <- 5251315 # RNG seed

# 2. Get data ----
source(here::here("analysis", "shelf_slope", "01_get_ss_data.R"))

cpue_comparison_df <- data.frame()

for(jj in 1:length(species_codes)) {
  
  sel_species <- species_codes[jj]
  
  dir.create(here::here("analysis", "shelf_slope", "output", sel_species), showWarnings = FALSE)
  
  # Setup effort data
  haul_df <- readRDS(file = here::here("analysis", "shelf_slope", "data", "ss_haul.rds"))
  
  area_swept_df <- dplyr::select(haul_df, MATCHUP, AREA_SWEPT_KM2, GEAR) |>
    dplyr::mutate(GEAR_COL = paste0("AREA_SWEPT_KM2_", GEAR)) |>
    dplyr::select(-GEAR) |>
    tidyr::pivot_wider(names_from = GEAR_COL, values_from = AREA_SWEPT_KM2)
  
  # Setup catch data
  catch_df <- readRDS(file = here::here("analysis", "shelf_slope", "data", "ss_catch.rds")) |>
    dplyr::filter(SPECIES_CODE == sel_species) |>
    dplyr::inner_join(dplyr::select(haul_df, HAULJOIN, MATCHUP, GEAR)) |>
    dplyr::select(WEIGHT, NUMBER_FISH, MATCHUP, GEAR) |>
    dplyr::mutate(NUMBER_FISH_COL = paste0("NUMBER_FISH_", GEAR),
                  WEIGHT_COL = paste0("WEIGHT_", GEAR)) |>
    dplyr:::select(-GEAR) |>
    dplyr::arrange(MATCHUP) |>
    tidyr::pivot_wider(names_from = c(NUMBER_FISH_COL, WEIGHT_COL),
                       values_from = c(NUMBER_FISH, WEIGHT),
                       values_fill = 0) |>
    dplyr::inner_join(area_swept_df)
  
  names(catch_df)[1:5] <- c("MATCHUP", "NUMBER_FISH_44", "NUMBER_FISH_172", "WEIGHT_44", "WEIGHT_172")
  
  catch_df <- catch_df |>
    dplyr::mutate(CPUE_WEIGHT_44 = WEIGHT_44/AREA_SWEPT_KM2_44,
                  CPUE_WEIGHT_172 = WEIGHT_172/AREA_SWEPT_KM2_172,
                  CPUE_NUMBER_FISH_44 = NUMBER_FISH_44/AREA_SWEPT_KM2_44,
                  CPUE_NUMBER_FISH_172 = NUMBER_FISH_172/AREA_SWEPT_KM2_172,
                  OFFSET = AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44,
                  RATIO = WEIGHT_44/WEIGHT_172)
  
  catch_mod_knots <- 10
  
  if(nrow(catch_df) < 10) {
    catch_mod_knots <- nrow(catch_df)-1
  }
  
  catch_mod <- gam(log(CPUE_WEIGHT_172+1) ~ s(log(CPUE_WEIGHT_44+1), 
                                              bs = "cr", 
                                              k = catch_mod_knots),
                   data = catch_df)
  
  catch_fit_df <- data.frame(CPUE_WEIGHT_44 = exp(seq(log(min(catch_df$CPUE_WEIGHT_44 + 1)),
                                                      log(max(catch_df$CPUE_WEIGHT_44 + 1)),
                                                      length = 100)) - 1)
  
  catch_fit_df$FIT_CPUE_WEIGHT_172 <- predict(catch_mod, 
                                             newdata = catch_fit_df)
  
  catch_fit_df$FIT_SE_CPUE_WEIGHT_172 <- predict(catch_mod, 
                                                newdata = catch_fit_df, se.fit = TRUE)$se
  
  # bias_correction_factor <- summary(catch_mod)$sigma^2/2
  bias_correction_factor <- 0
  
  catch_fit_df$FIT_LOWER_CPUE_WEIGHT_172 <- (exp(catch_fit_df$FIT_CPUE_WEIGHT_172 - catch_fit_df$FIT_SE_CPUE_WEIGHT_172 + bias_correction_factor)-1)
  catch_fit_df$FIT_UPPER_CPUE_WEIGHT_172 <- (exp(catch_fit_df$FIT_CPUE_WEIGHT_172 + catch_fit_df$FIT_SE_CPUE_WEIGHT_172 + bias_correction_factor)-1)
  catch_fit_df$FIT_CPUE_WEIGHT_172 <- (exp(catch_fit_df$FIT_CPUE_WEIGHT_172 + bias_correction_factor)-1)
  
  
  ragg::agg_png(file = here::here("analysis", "shelf_slope", "plots", paste0(sel_species, "_gam_ss.png")), width = 70, height = 70, units = "mm", res = 600)
  print(
    ggplot() +
      geom_ribbon(data = catch_fit_df,
                  mapping = aes(x = CPUE_WEIGHT_44,
                                ymin = FIT_LOWER_CPUE_WEIGHT_172,
                                ymax = FIT_UPPER_CPUE_WEIGHT_172),
                  alpha = 0.3,
                  color = NA) +
      geom_path(data = catch_fit_df,
                mapping = aes(x = CPUE_WEIGHT_44, 
                              y = FIT_CPUE_WEIGHT_172),
                color = "blue",
                linewidth = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      geom_point(data = dplyr::mutate(haul_df,
                                      YEAR = floor(CRUISE/100)) |> 
                   dplyr::select(YEAR, MATCHUP, VESSEL) |>
                   dplyr::inner_join(catch_df),
                 mapping = aes(x = CPUE_WEIGHT_44,
                               # shape = factor(YEAR)
                               y = CPUE_WEIGHT_172,
                 )) +
      theme_bw() +
      # scale_shape(name = "Year", solid = FALSE) +
      scale_x_continuous(name = expression(CPUE[30]~(kg %.%km^-2))) +
      scale_y_continuous(name = expression(CPUE[15]~(kg %.%km^-2)))
  )
  dev.off()
  
  cpue_comparison_df <- data.frame(SPECIES_CODE = sel_species) |>
    dplyr::bind_cols(as.data.frame(calculate_performance_metrics(x = catch_df$CPUE_WEIGHT_172, y = catch_df$CPUE_WEIGHT_44))) |>
    dplyr::bind_rows(cpue_comparison_df)
  
  
  # Setup length data
  length_df <- readRDS(file = here::here("analysis", "shelf_slope", "data", "ss_fish_crab_size.rds")) |>
    dplyr::filter(SPECIES_CODE == sel_species)
  
  if(nrow(length_df) < 0) {
    next
  }
  
  if(sum(length_df$FREQUENCY) < 150) {
    next
  }
  
  if(sel_species %in% c(68580, 68560)) {
    length_df$LENGTH <- length_df$WIDTH
  }
  
  # if(sel_species %in% c(68580, 68560, 69322, 69323)) {
  #   length_df$FREQUENCY <- 1
  # }
  
  # Setup length bins
  len_min <- min(length_df$LENGTH)
  len_max <- max(length_df$LENGTH)
  len_min <- len_min - len_min %% bin_width[jj]
  len_max <- len_max + (bin_width[jj] - len_max %% bin_width[jj])
  len_breaks <- seq(len_min, len_max, bin_width[jj])
  len_mid <- (len_breaks + bin_width[jj]/2)[1:(length(len_breaks) - 1)]
  
  length_df <- length_df |>
    dplyr::mutate(LEN_MIDPOINT = as.numeric(as.character(cut(LENGTH, breaks = len_breaks, labels = len_mid)))) |>
    dplyr::select(-SPECIES_CODE, -LENGTH) |>
    dplyr::group_by(HAULJOIN, LEN_MIDPOINT) |>
    dplyr::summarise(FREQUENCY = sum(FREQUENCY)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = "LEN_MIDPOINT", values_from = "FREQUENCY", values_fill = 0)
  
  length_df <- tidyr::pivot_longer(length_df, cols = 2:ncol(length_df)) |>
    dplyr::rename(LEN_MIDPOINT = name,
                  FREQUENCY = value) |>
    dplyr::inner_join(dplyr::select(haul_df, HAULJOIN, GEAR, MATCHUP)) |>
    dplyr::mutate(GEAR_COL = paste0("N_", GEAR)) |>
    dplyr::select(-HAULJOIN, -GEAR) |>
    tidyr::pivot_wider(names_from = GEAR_COL, values_from = FREQUENCY) |>
    dplyr::mutate(LEN_MIDPOINT = as.numeric(LEN_MIDPOINT)) |>
    dplyr::arrange(MATCHUP, LEN_MIDPOINT)
  
  # Set knots based on number of length bins
  gam_knots <- (length(len_mid)-1)-5 #  (NEED TO CHANGE WHEN THERE ARE MORE DATA)
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(sel_species == 69322) {
    gam_knots <- 5
  }
  
  pratio_df <- data.frame()
  unique_matchups <- unique(haul_df$MATCHUP)
  
  for(ii in 1:length(unique_matchups)) {
    
    sel_area_swept <- dplyr::filter(area_swept_df, MATCHUP == unique_matchups[ii])
    sel_catch <- dplyr::filter(catch_df, MATCHUP == unique_matchups[ii])
    sel_length <- dplyr::filter(length_df, MATCHUP == unique_matchups[ii])
    
    sel_length$p <- suppressMessages(
      selectivity_ratio(count1 = sel_length$N_44, 
                        count2 = sel_length$N_172, 
                        effort1 = sel_area_swept$AREA_SWEPT_KM2_44, 
                        effort2 = sel_area_swept$AREA_SWEPT_KM2_172)$p12
    )
    
    pratio_df <- dplyr::bind_rows(pratio_df, sel_length)
    
  }
  
  pratio_df <- dplyr::filter(pratio_df, !is.na(p))
  pratio_df$p_scaled <- sratio::scale_for_betareg(pratio_df$p, method = "sv")
  pratio_df$dummy_var <- 1
  
  
  # Setup four clusters and folds for each matchups
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  folds <- caret::groupKFold(group = pratio_df$MATCHUP)
  
  cv_results <- foreach::foreach(fold = folds) %dopar% {
    
    training_df <- pratio_df[fold, ]
    validation_df <- pratio_df[-fold, ]
    validation_df$dummy_var <- 0
    
    # Add in dummy station variable for predictions, to be added back in for output
    out_matchup <- validation_df$MATCHUP[1]
    validation_df$MATCHUP <- training_df$MATCHUP[1]
    
    gam_logit <- mgcv::gam(p_scaled ~ s(LEN_MIDPOINT, bs = "cr", k = gam_knots) + s(MATCHUP, bs = "re", by = dummy_var),
                           data = training_df |>
                             dplyr::mutate(MATCHUP = factor(MATCHUP)),
                           family = binomial(link = "logit"))
    
    gam_beta <- mgcv::gam(p_scaled ~ s(LEN_MIDPOINT, bs = "cr", k = gam_knots) + s(MATCHUP, bs = "re", by = dummy_var),
                          data = training_df |>
                            dplyr::mutate(MATCHUP = factor(MATCHUP)),
                          family = mgcv::betar(link = "logit"))
    
    fitted_logit <- predict(gam_logit, newdata = validation_df, type = "response")
    fitted_beta <- predict(gam_beta, newdata = validation_df, type = "response")
    
    validation_df$cv_fit_logit <- fitted_logit
    validation_df$cv_fit_beta <- fitted_beta
    
    # Reset matchup and dummy variable for fitting final models
    validation_df$MATCHUP <- out_matchup
    validation_df$dummy_var <- 1
    
    return(validation_df)
  }
  
  pratio_df <- do.call("rbind", cv_results)
  
  doParallel::stopImplicitCluster()
  
  # Fit models to generate mean prediction
  gam_logit <- gam(p_scaled ~ s(LEN_MIDPOINT, bs = "cr", k = gam_knots) + s(MATCHUP, bs = "re", by = dummy_var), 
                   data = pratio_df |>
                     dplyr::mutate(MATCHUP = factor(MATCHUP),
                                   dummy_var = 1),
                   family = binomial(link = "logit"))
  
  gam_beta <- gam(p_scaled ~ s(LEN_MIDPOINT, bs = "cr", k = gam_knots) + s(MATCHUP, bs = "re", by = dummy_var), 
                  data = pratio_df |>
                    dplyr::mutate(MATCHUP = factor(MATCHUP)),
                  family = betar(link = "logit"))
  
  logit_summary <- summary(gam_logit)
  beta_summary <- summary(gam_beta)
  
  rmse_df <- data.frame(model = c("Binomial", "Beta"), 
                        rmse = c((sum((pratio_df$cv_fit_logit-pratio_df$p_scaled)^2)/nrow(pratio_df))^0.5,
                                 (sum((pratio_df$cv_fit_beta-pratio_df$p_scaled)^2)/nrow(pratio_df))^0.5)
  )
  
  observed_prediction_df <- dplyr::bind_rows(
    make_prediction_df(model = gam_logit, lengths = len_min:len_max, model_name = "logit", type = "link"),
    make_prediction_df(model = gam_beta, lengths = len_min:len_max, model_name = "beta", type = "link")
  ) |>
    dplyr::mutate(type = "Observed")
  
  
  # Bootstrap to generate confidence intervals
  bootstrap_fits <- sratio_fit_bootstrap(dat = pratio_df, 
                                  lengths = len_min:len_max, 
                                  iterations = 1000, 
                                  clusters = 4, 
                                  seed = seed, 
                                  k = gam_knots)
  
  bootstrap_df <- do.call("rbind", bootstrap_fits)
  
  bootstrap_quantiles <- bootstrap_df |>
    dplyr::group_by(LEN_MIDPOINT) |>
    dplyr::summarise(q025 = quantile(fit, 0.025),
                     q975 = quantile(fit, 0.975),
                     q25 = quantile(fit, 0.25),
                     q75 = quantile(fit, 0.75)) |> 
    dplyr::mutate(p_q025 = inv_logit(q025),
                  sratio_q025 = 1/inv_logit(q025)-1,
                  p_q975 = inv_logit(q975),
                  sratio_q975 = 1/inv_logit(q975)-1,
                  p_q25 = inv_logit(q25),
                  sratio_q25 = 1/inv_logit(q25)-1,
                  p_q75 = inv_logit(q75),
                  sratio_q75 = 1/inv_logit(q75)-1) |>
    dplyr::mutate(type = "Observed")
  
  # Make plots of catch ratio and selectivity ratio ----
  plot_pratio <- ggplot() +
    geom_point(data = pratio_df,
               mapping = aes(x = LEN_MIDPOINT,
                             y = p_scaled),
               color = "grey70",
               shape = 21,
               size = rel(0.4)) +
    geom_ribbon(data = bootstrap_quantiles,
                mapping = aes(x = LEN_MIDPOINT,
                              ymin = p_q025,
                              max = p_q975),
                alpha = 0.5,
                fill = "grey20") +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = LEN_MIDPOINT,
                            y = p_q25),
              linetype = 3) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = LEN_MIDPOINT,
                            y = p_q75),
              linetype = 3) +
    geom_path(data = observed_prediction_df |>
                dplyr::filter(model == "logit"),
              mapping = aes(x = LEN_MIDPOINT,
                            y = p)) +
    scale_x_continuous(name = measurement_label[jj]) +
    scale_y_continuous(name = expression(italic(p['L,83-112,PNE']))) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  plot_obs_histogram <- ggplot() +
    geom_histogram(data = pratio_df,
                   mapping = aes(x = LEN_MIDPOINT),
                   bins = length(unique(pratio_df$LEN_MIDPOINT)),
                   fill = "grey70") +
    scale_x_continuous(name = measurement_label[jj], expand = c(0,0)) +
    scale_y_continuous(name = "Matchups (#)") +
    theme_bw()
  
  plot_sratio <- ggplot() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_ribbon(data = bootstrap_quantiles,
                mapping = aes(x = LEN_MIDPOINT,
                              ymin = sratio_q025,
                              max = sratio_q975),
                alpha = 0.5,
                fill = "grey20") +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = LEN_MIDPOINT,
                            y = sratio_q25),
              linetype = 3) +
    geom_path(data = bootstrap_quantiles,
              mapping = aes(x = LEN_MIDPOINT,
                            y = sratio_q75),
              linetype = 3) +
    geom_path(data = observed_prediction_df |>
                dplyr::filter(model == "logit"),
              mapping = aes(x = LEN_MIDPOINT,
                            y = sratio)) +
    scale_x_continuous(name = measurement_label[jj]) +
    scale_y_log10(name = expression(italic(S['L,PNE,83-112']))) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  # Write plots to file
  ragg::agg_png(file = here::here("analysis", "shelf_slope", "plots", paste0(sel_species, "_trawl_height_three_panel_ratios_n.png")), width = 169, height = 70, units = "mm", res = 600)
  print(cowplot::plot_grid(plot_pratio,
                           plot_obs_histogram,
                           plot_sratio,
                           nrow = 1,
                           labels = LETTERS[1:3]))
  dev.off()
  
  # Save species output to an rda file
  save(plot_pratio, 
       plot_obs_histogram, 
       plot_sratio, 
       rmse_df, 
       bootstrap_df, 
       bootstrap_fits, 
       observed_prediction_df, 
       haul_df, 
       length_df, 
       catch_df, 
       gam_logit, 
       gam_beta, 
       file = here::here("analysis", "shelf_slope", "output", sel_species, paste0("sratio_output_", sel_species, ".rda")))
}

write.csv(cpue_comparison_df,
          file = here::here("analysis", "shelf_slope", "plots", "cpue_comparison_metrics.csv"), 
          row.names = FALSE)
