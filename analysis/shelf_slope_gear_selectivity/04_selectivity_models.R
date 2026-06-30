library(ggthemes)
library(ggpp)
library(sratio)

source(here::here("analysis", "shelf_slope_gear_selectivity", "functions.R"))

# Load data ----
length_comp <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp.rds"))

length_comp_wide <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp_wide.rds"))

binned_comp <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "binned_length_comp.rds"))

binned_comp_wide <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "binned_length_comp_wide.rds"))

bootstrap_binned_comp_wide <- readRDS(here::here("analysis", "shelf_slope_gear_selectivity", "data", "bootstrap_binned_comp_wide.rds"))


analysis_species <- 
  data.frame(
    SPECIES_CODE = c(21740, 21720, 30060, 10110, 10112, 10130, 420, 435, 440, 455, 471, 472, 475, 477, 480, 485),
    COMMON_NAME = c(
      "walleye pollock", "Pacific cod", "Pacific ocean perch", "arrowtooth flounder", "Kamchatka flounder", "flathead sole",
      rep("skates", 10))
  )

gear_names <- 
  data.frame(GEAR = c(172, 44), GEAR_NAME = c("PNE", "83-112"))

# Length-frequency histograms ----

p_length_freq <- 
  ggplot() +
  geom_bar(
    data = dplyr::inner_join(length_comp, gear_names),
    mapping = aes(x = LENGTH, 
                  y = FREQUENCY*RAISING_FACTOR,
                  fill = GEAR_NAME),
    position = "dodge",
    stat = "identity",
    width = 1
  ) +
  ggpp::geom_text_npc(
    data = dplyr::select(length_comp, COMMON_NAME) |>
      unique(),
    mapping = aes(npcx = "right", npcy = "top", label = COMMON_NAME), 
    size = 3
  ) +
  scale_fill_colorblind(name = "Gear") +
  scale_y_continuous(name = "Frequency (#)", expand = expansion(mult = c(0,0.15)), limits = c(0, NA)) +
  scale_x_continuous(name = "Length (cm)", expand = c(0,0), limits = c(0, NA)) +
  facet_wrap(
    ~factor(COMMON_NAME, levels = unique(analysis_species$COMMON_NAME)), 
    scales = "free",
    ncol = 2
  ) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.77, 0.1),
        legend.direction = "horizontal",
        legend.title.position = "top",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "length_frequency.png"),
    units = "mm", res = 300, width = 120, height = 140)
print(p_length_freq)
dev.off()


p_binned_length_freq <- 
  ggplot() +
  geom_bar(
    data = dplyr::inner_join(binned_comp, gear_names),
    mapping = aes(x = LENGTH_BIN, 
                  y = FREQUENCY*RAISING_FACTOR,
                  fill = GEAR_NAME),
    position = "dodge",
    stat = "identity"
  ) +
  ggpp::geom_text_npc(
    data = dplyr::select(binned_comp, COMMON_NAME) |>
      unique(),
    mapping = aes(npcx = "right", npcy = "top", label = COMMON_NAME), 
    size = 3
  ) +
  scale_fill_colorblind(name = "Gear") +
  scale_y_continuous(name = "Frequency (#)", expand = expansion(mult = c(0,0.15)), limits = c(0, NA)) +
  scale_x_continuous(name = "Length (cm)", expand = c(0,0), limits = c(0, NA)) +
  facet_wrap(
    ~factor(COMMON_NAME, levels = unique(analysis_species$COMMON_NAME)), 
    scales = "free",
    ncol = 2
  ) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.77, 0.1),
        legend.direction = "horizontal",
        legend.title.position = "top",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "binned_length_frequency.png"),
    units = "mm", res = 300, width = 120, height = 140)
print(p_binned_length_freq)
dev.off()


# Fit selectivity models

unique_species <- unique(analysis_species$COMMON_NAME)

best_selectivity_model <- 
  selectivity_glmm_fits <- 
  model_tables <-
  bootstrap_fits <- 
  bootstrap_est <-
  species_plots <- 
  vector(mode = "list", length = length(unique_species))

prop_dat <- 
  binned_comp_wide |>
  dplyr::mutate(
    TOTAL = FREQUENCY_172 + FREQUENCY_44,
    P = FREQUENCY_172/(FREQUENCY_44+FREQUENCY_172),
    P_PLOT = FREQUENCY_172*RAISING_FACTOR_172/AREA_SWEPT_KM2_172/
      (FREQUENCY_172*RAISING_FACTOR_172/AREA_SWEPT_KM2_172+FREQUENCY_44*RAISING_FACTOR_44/AREA_SWEPT_KM2_44)
  )


for(ii in 1:length(unique_species)) {
  
  sel_spp <- unique_species[ii]
  
  sel_prop_dat <- dplyr::filter(prop_dat, COMMON_NAME == sel_spp)
  
  # Thygesen - Unbinned
  
  sel_dat <- length_comp |>
    dplyr::inner_join(
      dplyr::select(sratio::data_ss$haul, GEAR, MATCHUP, AREA_SWEPT_KM2)
    ) |>
    dplyr::filter(COMMON_NAME == sel_spp)

  length_values <- data.frame(MATCHUP = NA, LENGTH = seq(min(sel_dat$LENGTH), max(sel_dat$LENGTH), 4))

  sel_dat <- dplyr::bind_rows(sel_dat,
                              length_values) |>
    dplyr::arrange(LENGTH) |>
    tidyr::pivot_wider(names_from = "LENGTH", values_from = "FREQUENCY", values_fill = 0) |>
    dplyr::filter(!is.na(COMMON_NAME)) |>
    dplyr::arrange(GEAR, MATCHUP)

  input_unbinned <-
    list(
      N = as.matrix(sel_dat[, 6:ncol(sel_dat)]),
      SweptArea = sel_dat$AREA_SWEPT_KM2/sel_dat$RAISING_FACTOR,
      group = factor(sel_dat$MATCHUP),
      Gear = factor(sel_dat$GEAR, levels = c(172, 44)),
      Lvec = as.numeric(names(sel_dat)[6:ncol(sel_dat)])
    )

  lgcp_unbinned_mod1 <-
    sratio::gearcalib_fit(
      d = input_unbinned,
      model = "poisson"
    )

  # Thygesen - Binned
  sel_dat <- binned_comp |>
    dplyr::inner_join(
      dplyr::select(sratio::data_ss$haul, GEAR, MATCHUP, AREA_SWEPT_KM2)
    ) |>
    dplyr::filter(COMMON_NAME == sel_spp)
  # 
  # 
  # # Include dummy lengths
  # length_values <- data.frame(MATCHUP = NA, LENGTH_BIN = seq(min(sel_dat$LENGTH_BIN), max(sel_dat$LENGTH_BIN), 4))
  # 
  # sel_dat <- dplyr::bind_rows(sel_dat,
  #                             length_values) |>
  #   dplyr::arrange(LENGTH_BIN) |>
  #   tidyr::pivot_wider(names_from = "LENGTH_BIN", values_from = "FREQUENCY", values_fill = 0) |>
  #   dplyr::filter(!is.na(COMMON_NAME)) |>
  #   dplyr::arrange(GEAR, MATCHUP)
  # 
  # input_binned <-
  #   list(
  #     N = as.matrix(sel_dat[, 6:ncol(sel_dat)]),
  #     SweptArea = sel_dat$AREA_SWEPT_KM2/sel_dat$RAISING_FACTOR,
  #     group = factor(sel_dat$MATCHUP),
  #     Gear = factor(sel_dat$GEAR, levels = c(172, 44)),
  #     Lvec = as.numeric(names(sel_dat)[6:ncol(sel_dat)])
  #   )
  # 
  # lgcp_binned_mod1 <-
  #   sratio::gearcalib_fit(
  #     d = input_binned,
  #     model = "poisson"
  #   )
  # 
  # 
  # gearcalib_plot(lgcp_binned_mod1)
  
  # Binomial and betabinomial GAMMs
  
  
  bin_mod1 <- glmmTMB::glmmTMB(
    formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))),
    data = sel_prop_dat,
    weights = TOTAL,
    family = binomial(),
    verbose = FALSE
  )

  bin_mod2 <- glmmTMB::glmmTMB(
    formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))),
    data = sel_prop_dat,
    weights = TOTAL,
    family = binomial(),
    verbose = FALSE
  )
  
  bin_mod3 <- glmmTMB::glmmTMB(
    formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))),
    data = sel_prop_dat,
    weights = TOTAL,
    family = binomial(),
    verbose = FALSE
  )
  
  bb_mod1 <- glmmTMB::glmmTMB(
    formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod2 <- glmmTMB::glmmTMB(
    formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod3 <- glmmTMB::glmmTMB(
    formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod4 <- glmmTMB::glmmTMB(
    formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    dispformula = ~ LENGTH_BIN,
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod5 <- glmmTMB::glmmTMB(
    formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    dispformula = ~ LENGTH_BIN,
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod6 <- glmmTMB::glmmTMB(
    formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    dispformula = ~ LENGTH_BIN,
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod7 <- glmmTMB::glmmTMB(
    formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    dispformula = ~ s(LENGTH_BIN, bs = "tp", k = 5),
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod8 <- glmmTMB::glmmTMB(
    formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    dispformula = ~ s(LENGTH_BIN, bs = "tp", k = 5),
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  bb_mod9 <- glmmTMB::glmmTMB(
    formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_172 / RAISING_FACTOR_44))), 
    dispformula = ~ s(LENGTH_BIN, bs = "tp", k = 5),
    data = sel_prop_dat,
    weights = TOTAL,
    family = glmmTMB::betabinomial(),
    verbose = FALSE
  )
  
  selectivity_gamm_list <- 
    list(
      bin_mod1 = bin_mod1,
      bin_mod2 = bin_mod2,
      bin_mod3 = bin_mod3,
      bb_mod1 = bb_mod1,
      bb_mod2 = bb_mod2,
      bb_mod3 = bb_mod3,
      bb_mod4 = bb_mod4,
      bb_mod5 = bb_mod5,
      bb_mod6 = bb_mod6,
      bb_mod7 = bb_mod7,
      bb_mod8 = bb_mod8,
      bb_mod9 = bb_mod9
    )
  
  aic_table <- make_aic_table(selectivity_gamm_list)
  aic_table$COMMON_NAME <- sel_spp
  
  model_tables[[ii]] <- aic_table
  
  best_selectivity_model[[ii]] <- selectivity_gamm_list[[aic_table$model_name[1]]]
  
  plot(DHARMa::simulateResiduals(best_selectivity_model[[ii]]))
  
  best_fit <- init_fit <- data.frame(
    COMMON_NAME = sel_spp,
    LENGTH_BIN = seq(min(sel_prop_dat$LENGTH_BIN)-2, max(sel_prop_dat$LENGTH_BIN)+2, 1), 
    MODEL_NAME = aic_table$model_name[1],
    AREA_SWEPT_KM2_44 = 1, 
    AREA_SWEPT_KM2_172 = 1, 
    RAISING_FACTOR_44 = 1, 
    RAISING_FACTOR_172 = 1, 
    TOTAL = 1,
    MATCHUP = factor(-999))
  
  best_fit[, c('link_fit', 'link_se')] <- 
    predict(best_selectivity_model[[ii]], 
            type = "link", 
            newdata = best_fit,
            se.fit = TRUE,
            re.form = NA,
            allow.new.levels = TRUE)
  
  best_fit$fit <- best_selectivity_model[[ii]]$modelInfo$family$linkinv(best_fit$link_fit)
  best_fit$ci_lwr <- best_selectivity_model[[ii]]$modelInfo$family$linkinv(best_fit$link_fit - 2*best_fit$link_se)
  best_fit$ci_upr <- best_selectivity_model[[ii]]$modelInfo$family$linkinv(best_fit$link_fit + 2*best_fit$link_se) 
  
  selectivity_glmm_fits[[ii]] <- best_fit
  
  # Bootstrap fits for most parsimonious model ----
  sel_bootstrap <- bootstrap_binned_comp_wide[[sel_spp]]
  
  start_time <- Sys.time()
  message(start_time, ": Start bootstrap fits")
  bootstrap_fits[[ii]] <- lapply(
    X = sel_bootstrap,
    FUN = function(x, mod, fit_df) {
      x_sel <- x |>
        dplyr::mutate(
          TOTAL = FREQUENCY_172 + FREQUENCY_44,
          P = FREQUENCY_172/(FREQUENCY_44+FREQUENCY_172)
        )
      
      x_mod <- update(mod, data = x_sel)
      
      if(x_mod$fit$convergence == 0 & x_mod$sdr$pdHess == TRUE) {
        fit_df[, c('link_fit', 'link_se')] <- 
          predict(x_mod, 
                  type = "link", 
                  newdata = fit_df,
                  se.fit = TRUE,
                  re.form = NA,
                  allow.new.levels = TRUE)
        
        fit_df$fit <- x_mod$modelInfo$family$linkinv(fit_df$link_fit)
        # fit_df$ci_lwr <- x_mod$modelInfo$family$linkinv(fit_df$link_fit - 2*fit_df$link_se)
        # fit_df$ci_upr <- x_mod$modelInfo$family$linkinv(fit_df$link_fit + 2*fit_df$link_se) 
        
        return(fit_df)
      } else {
        return(NULL)
      }

    },
    mod = best_selectivity_model[[ii]],
    fit_df = init_fit
    
    )
  
  message(Sys.time(), ": End bootstrap fits (", round(difftime(Sys.time(), start_time, units = "mins"), 1), " min)")
  
  bootstrap_est[[ii]] <- 
    do.call(what = dplyr::bind_rows, bootstrap_fits[[ii]]) |>
    dplyr::group_by(COMMON_NAME, LENGTH_BIN) |>
    dplyr::summarise(P_CI_LWR = quantile(fit, 0.025),
                     P_CI_UPR = quantile(fit, 0.975),
                     P_MEDIAN = median(fit),
                     P_25 = quantile(fit, 0.25),
                     P_75 = quantile(fit, 0.75))
  
  # Plot the most parsimonious model (no bootstrap)
  p_species_plot <- 
    ggplot() +
    geom_point(
      data = sel_prop_dat, 
      mapping = aes(
        x = LENGTH_BIN, 
        y = P_PLOT, 
        size = cut(TOTAL, c(0, 10, 20, 50, Inf))
      ), 
      alpha = 0.4, 
      color = "tan") +
    geom_ribbon(data = best_fit, mapping = aes(x = LENGTH_BIN, ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
    geom_path(data = best_fit, mapping = aes(x = LENGTH_BIN, y = fit)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "purple4") +
    scale_y_continuous(name = "Relative selectivity", limits = c(0, 1)) +
    scale_x_continuous(name = "Length (cm)", expand = c(0,0)) +
    scale_size_discrete(name = "# Lengths") +
    theme_bw()
  
  # Plot the most parsimonious model (bootstrap)
  p_species_boot_plot <- 
    ggplot() +
    geom_point(
      data = sel_prop_dat, 
      mapping = aes(
        x = LENGTH_BIN, 
        y = P_PLOT, 
        size = cut(TOTAL, c(0, 10, 20, 50, Inf))
      ), 
      alpha = 0.4, 
      color = "tan") +
    geom_ribbon(data = bootstrap_est[[ii]], mapping = aes(x = LENGTH_BIN, ymin = P_CI_LWR, ymax = P_CI_UPR), alpha = 0.3) +
    geom_path(data = bootstrap_est[[ii]], mapping = aes(x = LENGTH_BIN, y = P_MEDIAN)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "purple4") +
    scale_y_continuous(name = "Relative selectivity", limits = c(0, 1)) +
    scale_x_continuous(name = "Length (cm)", expand = c(0,0)) +
    scale_size_discrete(name = "# Lengths") +
    theme_bw()
  
  species_plots[[ii]] <-
  cowplot::plot_grid(
    p_species_plot + ggtitle(paste0(sel_spp, " raw fit")), 
                     p_species_boot_plot + ggtitle(paste0(sel_spp, " bootstrap")),
                     ncol = 2,
                     align = "hv"
    )
  gearcalib_plot(lgcp_unbinned_mod1)
  print(species_plots[[ii]])
  
}

saveRDS(
  object = species_plots, 
  here::here("analysis", "shelf_slope_gear_selectivity", "output", "species_plots.rds")
)
saveRDS(
  object = bootstrap_est, 
  here::here("analysis", "shelf_slope_gear_selectivity", "output", "bootstrap_est.rds")
)
saveRDS(
  object = bootstrap_fits, 
  here::here("analysis", "shelf_slope_gear_selectivity", "output", "bootstrap_fits.rds")
)
saveRDS(
  object = selectivity_glmm_fits,
  here::here("analysis", "shelf_slope_gear_selectivity", "output", "selectivity_glmm_fits.rds")
)
saveRDS(
  object = model_tables,
  here::here("analysis", "shelf_slope_gear_selectivity", "output", "selectivity_model_tables.rds")
)



bootstrap_plot_data <- do.call(what = dplyr::bind_rows, bootstrap_est)

(p_glmm_selectivity <- 
    ggplot() +
    geom_point(
      data = prop_dat, 
      mapping = 
        aes(
          x = LENGTH_BIN, 
          y = P_PLOT, 
          size = cut(TOTAL, c(0, 10, 20, 50, Inf))
        ), 
      alpha = 0.2, 
      color = "tan"
    ) +
    geom_ribbon(
      data = bootstrap_plot_data, 
      mapping = 
        aes(
          x = LENGTH_BIN, 
          ymin = P_CI_LWR, 
          ymax = P_CI_UPR
        ), 
      alpha = 0.3
    ) +
    geom_path(data = bootstrap_plot_data, mapping = aes(x = LENGTH_BIN, y = P_MEDIAN)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "purple4") +
    scale_y_continuous(name = "Relative catch efficiency", limits = c(0, 1)) +
    scale_x_continuous(name = "Length (cm)", expand = c(0,0)) +
    scale_size_discrete(name = "# Lengths") +
    theme_bw() +
    facet_wrap(~factor(COMMON_NAME, levels = unique_species), scales = "free_x", ncol = 2))

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "selectivity_ratio_glmm.png"),
    units = "mm", res = 300, width = 169, height = 140)
print(p_glmm_selectivity)
dev.off()


