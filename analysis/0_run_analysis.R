# Selectivity ratio for 15/30 comparison

# 1. Setup ----
library(sratio)

n_cores <- 4 # Cores for parallel processing
bin_width <- 5 # May need to change for different species
species_codes <- c(21740, 21720, 10210, 10261, 10110, 10112, 10130, 10285, 471, 68580, 68560, 69322, 69323)
seed <- 909823 # RNG seed

# 2. Get data ----
source(here::here("analysis", "1_get_data.R"))


for(jj in 1:length(species_codes)) {
  
  sel_species <- species_codes[jj]
  
  dir.create(here::here("output", sel_species), showWarnings = FALSE)
  
  # Setup area spet data
  haul_df <- readRDS(file = here::here("data", "hauls_1530.rds"))
  
  area_swept_df <- dplyr::select(haul_df, MATCHUP, AREA_SWEPT_KM2, TREATMENT) |>
    dplyr::mutate(TREATMENT_COL = paste0("AREA_SWEPT_KM2_", TREATMENT)) |>
    dplyr::select(-TREATMENT) |>
    tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = AREA_SWEPT_KM2)
  
  # Setup catch data
  catch_df <- readRDS(file = here::here("data", "catch_1530.rds")) |>
    dplyr::filter(SPECIES_CODE == sel_species) |>
    dplyr::inner_join(dplyr::select(haul_df, HAULJOIN, MATCHUP))
  
  
  # Setup length data
  length_df <- readRDS(file = here::here("data", "fish_lengths_1530.rds")) |>
    dplyr::filter(SPECIES_CODE == sel_species)
  
  # Setup length bins
  len_min <- min(length_df$LENGTH)
  len_max <- max(length_df$LENGTH)
  len_min <- len_min - len_min %% bin_width
  len_max <- len_max + (bin_width - len_max %% bin_width)
  len_breaks <- seq(len_min, len_max, bin_width)
  len_mid <- (len_breaks + bin_width/2)[1:(length(len_breaks) - 1)]
  
  gam_knots <- (length(len_mid)-1)
  
  if(gam_knots > 10) {
    gam_knots <- 10
  }
  
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
    dplyr::inner_join(dplyr::select(haul_df, HAULJOIN, TREATMENT, MATCHUP)) |>
    dplyr::mutate(TREATMENT_COL = paste0("N_", TREATMENT)) |>
    dplyr::select(-HAULJOIN, -TREATMENT) |>
    tidyr::pivot_wider(names_from = TREATMENT_COL, values_from = FREQUENCY) |>
    dplyr::mutate(LEN_MIDPOINT = as.numeric(LEN_MIDPOINT)) |>
    dplyr::arrange(MATCHUP, LEN_MIDPOINT)
  
  pratio_df <- data.frame()
  unique_matchups <- unique(haul_df$MATCHUP)
  
  for(ii in 1:length(unique_matchups)) {
    
    sel_area_swept <- dplyr::filter(area_swept_df, MATCHUP == unique_matchups[ii])
    sel_catch <- dplyr::filter(catch_df, MATCHUP == unique_matchups[ii])
    sel_length <- dplyr::filter(length_df, MATCHUP == unique_matchups[ii])
    
    sel_length$p <- selectivity_ratio(count1 = sel_length$N_30, 
                                      count2 = sel_length$N_15, 
                                      effort1 = sel_area_swept$AREA_SWEPT_KM2_30, 
                                      effort2 = sel_area_swept$AREA_SWEPT_KM2_15)$p12
    
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
  bootstrap_fits <- run_bootstrap(dat = pratio_df, 
                                  lengths = len_min:len_max, 
                                  iterations = 1000, 
                                  clusters = 4, 
                                  seed = seed, 
                                  k = gam_knots)
  
  bootstrap_df <- do.call("rbind", bootstrap_fits)
  
  bootstrap_quantiles <- bootstrap_df |>
    dplyr::group_by(LEN_MIDPOINT) |>
    dplyr::summarise(q025 = quantile(fit, 0.025),
                     q975 = quantile(fit, 0.975)) |> 
    dplyr::mutate(p_q025 = inv_logit(q025),
                  sratio_q025 = 1/inv_logit(q025)-1,
                  p_q975 = inv_logit(q975),
                  sratio_q975 = 1/inv_logit(q975)-1) |>
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
    geom_path(data = observed_prediction_df |>
                dplyr::filter(model == "logit"),
              mapping = aes(x = LEN_MIDPOINT,
                            y = p)) +
    scale_x_continuous(name = "Fork Length (cm)") +
    scale_y_continuous(name = expression(italic(p['L,30,15']))) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  plot_obs_histogram <- ggplot() +
    geom_histogram(data = pratio_df,
                   mapping = aes(x = LEN_MIDPOINT),
                   bins = length(unique(pratio_df$LEN_MIDPOINT)),
                   fill = "grey70") +
    scale_x_continuous(name = "Fork Length (cm)", expand = c(0,0)) +
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
    geom_path(data = observed_prediction_df |>
                dplyr::filter(model == "logit"),
              mapping = aes(x = LEN_MIDPOINT,
                            y = sratio)) +
    scale_x_continuous(name = "Fork Length (cm)") +
    scale_y_log10(name = expression(italic(S['L,15,30']))) +
    scale_color_tableau() +
    scale_fill_tableau() +
    theme_bw()
  
  # Write plots to file
  png(file = here::here("plots", paste0(sel_species, "_trawl_height_three_panel_ratios_n.png")), width = 169, height = 70, units = "mm", res = 600)
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
       file = here::here("output", sel_species, paste0("sratio_output_", sel_species, ".rda")))
}
