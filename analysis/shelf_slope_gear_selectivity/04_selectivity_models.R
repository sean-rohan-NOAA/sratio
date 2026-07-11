library(ggthemes)
library(ggpp)
library(sratio)
library(glmmTMB)

source(here::here("analysis", "shelf_slope_gear_selectivity", "functions.R"))

# Load data ----
length_comp <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp.rds"))

length_comp_wide <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "length_comp_wide.rds"))

binned_comp <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "binned_length_comp.rds"))

binned_comp_wide <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "binned_length_comp_wide.rds"))

prop_dat <- 
  binned_comp_wide |>
  dplyr::mutate(
    TOTAL = FREQUENCY_172 + FREQUENCY_44,
    P = FREQUENCY_172/(FREQUENCY_44+FREQUENCY_172),
    P_PLOT = FREQUENCY_172*RAISING_FACTOR_172/AREA_SWEPT_KM2_172/
      (FREQUENCY_172*RAISING_FACTOR_172/AREA_SWEPT_KM2_172+FREQUENCY_44*RAISING_FACTOR_44/AREA_SWEPT_KM2_44)
  )

catch_data <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "data", "catch_data.rds"))

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

unique_species <- unique(analysis_species$COMMON_NAME)


# Bootstrap aggregate relative catch efficiency ratios

aggregate_catch_eff <- vector(mode = "list", length = length(bootstrap_binned_comp_wide))

for(hh in 1:length(bootstrap_binned_comp_wide)) {
  
    catch_eff <- 
    lapply(X = bootstrap_binned_comp_wide[[hh]],
           FUN = function(x) {
             x |>
               dplyr::mutate(
                 N_44 = FREQUENCY_44*RAISING_FACTOR_44,
                 N_172 = FREQUENCY_172*RAISING_FACTOR_172
               ) |>
               dplyr::group_by(LENGTH_BIN, COMMON_NAME) |>
               dplyr::summarise(
                 TOTAL_N_44 = sum(N_44),
                 TOTAL_N_172 = sum(N_172),
                 TOTAL_AREA_SWEPT_KM2_44 = sum(AREA_SWEPT_KM2_44),
                 TOTAL_AREA_SWEPT_KM2_172 = sum(AREA_SWEPT_KM2_44),
                 .groups = "keep"
               ) |>
               dplyr::ungroup() |>
               dplyr::mutate(
                 TOTAL_CPUE_44 = TOTAL_N_44/TOTAL_AREA_SWEPT_KM2_44,
                 TOTAL_CPUE_172 = TOTAL_N_172/TOTAL_AREA_SWEPT_KM2_172,
                 BOOT_P = TOTAL_CPUE_172/(TOTAL_CPUE_44+TOTAL_CPUE_172))
           }) 
  
  aggregate_catch_eff[[hh]] <- 
    do.call(what = dplyr::bind_rows, catch_eff) |>
    dplyr::group_by(COMMON_NAME, LENGTH_BIN) |>
    dplyr::summarise(
      BOOT_P_MEAN = mean(BOOT_P),
      BOOT_P_MEDIAN = median(BOOT_P),
      BOOT_P_CI_LWR = quantile(BOOT_P, 0.025),
      BOOT_P_CI_UPR = quantile(BOOT_P, 0.975),
      BOOT_P_Q25 = quantile(BOOT_P, 0.25),
      BOOT_P_Q75 = quantile(BOOT_P, 0.75),
      .groups = "keep"
    ) |>
    dplyr::ungroup()
  
}

aggregate_catch_eff <-  do.call(what = dplyr::bind_rows, aggregate_catch_eff)

saveRDS(aggregate_catch_eff, file = here::here("analysis", "shelf_slope_gear_selectivity", "output", "aggregate_catch_eff.rds"))


# Catch by haul pair ----
p_numerical_catch_by_pair <- 
  ggplot() +
  geom_bar(data = dplyr::inner_join(catch_data, gear_names) |>
             dplyr::inner_join(dplyr::select(sratio::data_ss$haul, HAULJOIN, AREA_SWEPT_KM2)) |>
             dplyr::arrange(-GEAR),
           mapping = aes(x = MATCHUP, y = NUMBER_FISH, fill = GEAR_NAME, group = MATCHUP),
           stat = "identity", position = "stack") +
  ggpp::geom_text_npc(
    data = dplyr::select(length_comp, COMMON_NAME) |>
      unique(),
    mapping = aes(npcx = "left", npcy = "top", label = COMMON_NAME), 
    size = 3
  ) +
  scale_y_continuous(name = "Frequency", expand = expansion(mult = c(0,0.15)), limits = c(0, NA)) +
  scale_x_continuous(name = "Tow pair") +
  scale_fill_colorblind(name = "Gear") +
  facet_wrap(~factor(COMMON_NAME, levels = unique_species), scales = "free_y", ncol = 2) +
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


p_numerical_cpue_by_pair <- 
  ggplot() +
  geom_point(data = dplyr::inner_join(catch_data, gear_names) |>
             dplyr::inner_join(dplyr::select(sratio::data_ss$haul, HAULJOIN, AREA_SWEPT_KM2)),
           mapping = aes(x = MATCHUP, y = NUMBER_FISH/AREA_SWEPT_KM2, color = GEAR_NAME, group = MATCHUP)) +
  ggpp::geom_text_npc(
    data = dplyr::select(length_comp, COMMON_NAME) |>
      unique(),
    mapping = aes(npcx = "left", npcy = "top", label = COMMON_NAME), 
    size = 3
  ) +
  scale_y_log10(name = expression(CPUE*' (#/km'^2*')'), expand = c(0, NA)) +
  scale_x_continuous(name = "Tow pair") +
  scale_color_colorblind(name = "Gear") +
  facet_wrap(~factor(COMMON_NAME, levels = unique_species), scales = "free_y", ncol = 2) +
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

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "numerical_cpue_by_pair.png"),
    units = "mm", res = 300, width = 120, height = 140)
print(p_numerical_cpue_by_pair)
dev.off()

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "numerical_catch_by_pair.png"),
    units = "mm", res = 300, width = 120, height = 140)
print(p_numerical_catch_by_pair)
dev.off()

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

best_selectivity_model <- 
  selectivity_glmm_fits <- 
  model_tables <-
  bootstrap_fits <- 
  bootstrap_est <-
  species_plots <- 
  vector(mode = "list", length = length(unique_species))




for(ii in 1:length(unique_species)) {
  
  sel_spp <- unique_species[ii]
  
  sel_prop_dat <- dplyr::filter(prop_dat, COMMON_NAME == sel_spp)
  
  # Binomial and betabinomial GAMMs
  
  bin_mod1 <- 
    try(
      glmmTMB::glmmTMB(
        formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))),
        data = sel_prop_dat,
        weights = TOTAL,
        family = binomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bin_mod2 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))),
        data = sel_prop_dat,
        weights = TOTAL,
        family = binomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bin_mod3 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))),
        data = sel_prop_dat,
        weights = TOTAL,
        family = binomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod1 <-
    try(
      glmmTMB::glmmTMB(
        formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod2 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod3 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod4 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        dispformula = ~ LENGTH_BIN,
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod5 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        dispformula = ~ LENGTH_BIN,
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod6 <- 
    try(
      glmmTMB::glmmTMB(
        formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        dispformula = ~ LENGTH_BIN,
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod7 <- 
    try( 
      glmmTMB::glmmTMB(
        formula = P ~ 1 + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        dispformula = ~ s(LENGTH_BIN, bs = "tp", k = 5),
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod8 <- 
    try(
      glmmTMB::glmmTMB(
        formula = P ~ LENGTH_BIN + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        dispformula = ~ s(LENGTH_BIN, bs = "tp", k = 5),
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  bb_mod9 <-  
    try(
      glmmTMB::glmmTMB(
        formula = P ~ s(LENGTH_BIN, bs = "tp", k = 5) + (1 | MATCHUP) + offset(I(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44  * RAISING_FACTOR_44 / RAISING_FACTOR_172))), 
        dispformula = ~ s(LENGTH_BIN, bs = "tp", k = 5),
        data = sel_prop_dat,
        weights = TOTAL,
        family = glmmTMB::betabinomial(),
        verbose = FALSE
      ), 
      silent = TRUE)
  
  selectivity_model_list <- 
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
  
  # Remove models that could not be fit - not many
  selectivity_model_list <- selectivity_model_list[lapply(selectivity_model_list, FUN = is, class2 = "try-error") == 0]
  
  aic_table <- make_aic_table(selectivity_model_list)
  aic_table$COMMON_NAME <- sel_spp
  
  model_tables[[ii]] <- aic_table
  
  best_model_name <- aic_table$model_name[aic_table$best_model]

  best_selectivity_model[[ii]] <- selectivity_model_list[[best_model_name]]
  
  # plot(DHARMa::simulateResiduals(best_selectivity_model[[ii]]))
  
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
      
      x_mod <- try(update(mod, data = x_sel), silent = TRUE)
      
      if(is(x_mod, "try-error")) {
        return(NULL)
      }
      
      if(x_mod$fit$convergence == 0 & x_mod$sdr$pdHess == TRUE) {
        fit_df[, c('link_fit', 'link_se')] <- 
          predict(x_mod, 
                  type = "link", 
                  newdata = fit_df,
                  se.fit = TRUE,
                  re.form = NA,
                  allow.new.levels = TRUE)
        
        fit_df$fit <- x_mod$modelInfo$family$linkinv(fit_df$link_fit)
        
        return(fit_df)
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
  # gearcalib_plot(lgcp_unbinned_mod1)
  print(species_plots[[ii]])
  
}

# Make model tables

# top_models  <- lapply(X = model_tables,
#                       FUN = function(x) {
#                         x_pass_check <- x[x$pass_check == TRUE, ]
#                         x_pass_check[1:3,]
#                       })
# 
# top_models <- do.call(what = dplyr::bind_rows, args = top_models) |>
#   dplyr::select(COMMON_NAME, model_name, aic, k, convergence, max_gradient, delta_aic)
# 
# write.csv(
#   top_models,
#   here::here("analysis", "shelf_slope_gear_selectivity", "plots", "top_selectivity_models.csv"),
#   row.names = FALSE
# )


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


# Plot selectivity results

species_plots <- readRDS(here::here("analysis", "shelf_slope_gear_selectivity", "output", "species_plots.rds"))
bootstrap_est <- readRDS(here::here("analysis", "shelf_slope_gear_selectivity", "output", "bootstrap_est.rds"))
bootstrap_fits <- readRDS(here::here("analysis", "shelf_slope_gear_selectivity", "output", "bootstrap_fits.rds"))
selectivity_glmm_fits <- readRDS(here::here("analysis", "shelf_slope_gear_selectivity", "output", "selectivity_glmm_fits.rds"))
model_tables <- readRDS(here::here("analysis", "shelf_slope_gear_selectivity", "output", "selectivity_model_tables.rds"))
aggregate_catch_eff <- readRDS(file = here::here("analysis", "shelf_slope_gear_selectivity", "output", "aggregate_catch_eff.rds"))

bootstrap_plot_data <- do.call(what = dplyr::bind_rows, bootstrap_est)

p_glmm_selectivity <- 
    ggplot() +
    geom_point(
      data = prop_dat, 
      mapping = 
        aes(
          x = LENGTH_BIN, 
          y = P_PLOT, 
          size = cut(TOTAL, c(0, 10, 20, 50, Inf), labels = c("1-10", "11-20", "21-50", ">50"))
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
      alpha = 0.3, fill = "blue"
    ) +
    geom_path(data = bootstrap_plot_data, mapping = aes(x = LENGTH_BIN, y = P_MEDIAN), color = "blue") +
    geom_point(data = aggregate_catch_eff,
              mapping = aes(x = LENGTH_BIN, y = BOOT_P_MEDIAN)) +
    geom_path(data = aggregate_catch_eff,
              mapping = aes(x = LENGTH_BIN, y = BOOT_P_CI_LWR), linetype = 3) +
    geom_path(data = aggregate_catch_eff,
              mapping = aes(x = LENGTH_BIN, y = BOOT_P_CI_UPR), linetype = 3) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "purple4") +
    facet_wrap(~factor(COMMON_NAME, levels = unique_species), scales = "free_x", ncol = 2) +
    scale_y_continuous(name = "Relative catch efficiency", limits = c(0, 1)) +
    scale_x_continuous(name = "Length (cm)", expand = c(0,0)) +
    scale_size_discrete(name = "Lengths") +
    scale_fill_manual(name = NULL, values = c("Best-fit model" = "blue", "Pop.-level bootstrap" = "black"), guide = FALSE) +
    scale_color_manual(name = NULL, values = c("Best-fit model" = "blue", "Pop.-level bootstrap" = "black"), guide = FALSE) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9, hjust = 0, face = "bold"),
          legend.position = "inside",
          legend.position.inside = c(0.77, 0.1),
          legend.direction = "horizontal",
          legend.title.position = "top",
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9))

png(here::here("analysis", "shelf_slope_gear_selectivity", "plots", "selectivity_ratio_glmm.png"),
    units = "mm", res = 300, width = 169, height = 140)
print(p_glmm_selectivity)
dev.off()

# ggplot() +
#   geom_point(
#     data = prop_dat, 
#     mapping = 
#       aes(
#         x = LENGTH_BIN, 
#         y = P_PLOT, 
#         size = cut(TOTAL, c(0, 10, 20, 50, Inf))
#       ), 
#     alpha = 0.2, 
#     color = "tan"
#   ) +
#   geom_ribbon(
#     data = bootstrap_plot_data, 
#     mapping = 
#       aes(
#         x = LENGTH_BIN, 
#         ymin = P_CI_LWR, 
#         ymax = P_CI_UPR
#       ), 
#     alpha = 0.3
#   ) +
#   geom_path(data = bootstrap_plot_data, mapping = aes(x = LENGTH_BIN, y = P_MEDIAN)) +
#   geom_hline(yintercept = 0.5, linetype = 2, color = "purple4") +
#   scale_y_continuous(name = "Relative catch efficiency", limits = c(0, 1)) +
#   scale_x_continuous(name = "Length (cm)", expand = c(0,0)) +
#   scale_size_discrete(name = "# Lengths") +
#   theme_bw() +
#   facet_wrap(~factor(COMMON_NAME, levels = unique_species), scales = "free_x", ncol = 2)

# Make selectivity model AIC table ----
full_aic_table <- 
  model_tables |>
  do.call(what = dplyr::bind_rows) |>
  dplyr::filter(pass_check == TRUE) |>
  dplyr::select(COMMON_NAME, model_name, delta_aic) |>
  dplyr::mutate(delta_aic = format(delta_aic, nsmall = 2),
                model_name = stringr::str_replace(string = model_name, pattern = "_mod", replacement = ""),
                model_name = stringr::str_replace(string = model_name, pattern = "bin", replacement = "bi"),
                model_name = toupper(model_name)
                )|>
  dplyr::arrange(model_name) |>
  tidyr::pivot_wider(values_from = "delta_aic", names_from = "model_name", values_fill = "")

  model_tables |>
  do.call(what = dplyr::bind_rows) |>
  dplyr::filter(best_model == TRUE) |>
    dplyr::select(COMMON_NAME, model_name, delta_aic)
  
  model_tables |>
    do.call(what = dplyr::bind_rows) |>
    dplyr::select(model_name, k) |>
    unique() |>
    dplyr::arrange(model_name)

write.csv(
  full_aic_table, 
  here::here("analysis", "shelf_slope_gear_selectivity", "plots", 
             "selectivity_aicc_full_table.csv"), 
  row.names = FALSE)
