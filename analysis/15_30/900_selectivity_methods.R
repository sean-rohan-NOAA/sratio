# Compare selectivity ratios using Thygesen and Kotwicki methods

library(sratio) # Implements Thygesen, Kotwicki, and other methods
library(selfisher)
library(splines)
library(cowplot)
library(scales)

# spp_code <- 10210
xlab <- unique(sratio::species_code_label(x = spp_code))
common_name <- unique(sratio::species_code_label(x = spp_code, type = "common_name"))

# Setup directories if they don't exist ------------------------------------------------------------

dir.create(
  here::here("analysis", "15_30", "plots", "selectivity_ratios"), 
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here("analysis", "15_30", "output", "selectivity_ratios"), 
  recursive = TRUE,
  showWarnings = FALSE
)

# Data ---- ----------------------------------------------------------------------------------------
catch_at_size_dat <- readRDS(file = here::here("analysis", "15_30", "output", "catch_at_length_1530.rds")) |>
  dplyr::filter(SPECIES_CODE %in% spp_code)

sratio_dat <- 
  readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
  dplyr::mutate(N_TOTAL = N_15+N_30,
                p12 = N_30/N_TOTAL) |>
  dplyr::filter(SPECIES_CODE %in% spp_code,
                N_TOTAL > 0)

boot_dat <- 
  readRDS(here::here("analysis", "15_30", "output", spp_code, paste0("bootstrap_samples_", spp_code, ".rds")))

# Number of bootstrap samples to use -- manually reduce when testing
n_boot <- length(boot_dat$wide)

# Themes and colors --------------------------------------------------------------------------------

year_colors <- c(`1995` = "#0072B2", 
                 `1998` =  "#F0E442", 
                 `2021` =  "#009E73", 
                 `2022` =  "#56B4E9", 
                 `2023` = "#000000", 
                 `2024` = "#E69F00")


# Thygesen et al. (2019) log-Gaussian Cox process model --------------------------------------------
lgcp_dat <-
  catch_at_size_dat |>
  dplyr::mutate(TOTAL_COUNT = round(FREQUENCY * SAMPLING_FACTOR)) |>
  dplyr::select(-SAMPLING_FACTOR, -FREQUENCY) |>
  dplyr::group_by(HAULJOIN, MATCHUP, TREATMENT, SPECIES_CODE, AREA_SWEPT_KM2, SIZE_BIN) |>
  dplyr::summarise(TOTAL_COUNT = sum(TOTAL_COUNT)) |>
  dplyr::arrange(SIZE_BIN) |>  
  tidyr::pivot_wider(values_from = "TOTAL_COUNT", names_from = "SIZE_BIN", values_fill = 0)

lgcp_input <- 
  list(
    N = as.matrix(lgcp_dat[, 6:ncol(lgcp_dat)]),
    SweptArea = lgcp_dat$AREA_SWEPT_KM2,
    group = factor(lgcp_dat$MATCHUP),
    Gear = factor(lgcp_dat$TREATMENT),
    Lvec = as.numeric(names(lgcp_dat)[6:ncol(lgcp_dat)])
  )

# Bootstrap estimate of aggregate CPUE ratio - use throughout
boot_results <- gearcalib_boot(lgcp_input, quantiles = c(0.025,0.5,0.975), nboot = 1000)

survey_level_bootstrap <- 
  boot_results$BootQuantiles |>
  t() |>
  data.frame() |>
  dplyr::mutate(
    s_mean = boot_results$RawEstimate,
    SIZE_BIN = as.numeric(colnames(boot_results$BootQuantiles))
    )

names(survey_level_bootstrap)[1:3] <-  c("s_q025", "s_q500", "s_q975")


# Fit LGCP model
use_logsd <- TRUE
lgcp_fit <- gearcalib_fit(d = lgcp_input, model = "poisson")

# Fix random walk logsd at a tiny value when the Hessian is not positive definite
if(!lgcp_fit$rep$pdHess) {
  lgcp_fit <- gearcalib_fit(d = lgcp_input, model = "poisson", logsdGearRW = -10)
  use_logsd <- FALSE
}

lgcp_plots <- 
  gearcalib_plot(
    fit = lgcp_fit, 
    boot = boot_results, 
    add_bootquantiles = TRUE, 
    xlab = xlab
  )

n_hauls <-
  catch_at_size_dat|>
  dplyr::filter(SPECIES_CODE %in% spp_code) |>
  dplyr::group_by(TREATMENT, SIZE_BIN) |>
  dplyr::summarise(n_positive = n(),
                   .groups = "keep")

p_encounters <- 
  ggplot() +
  geom_path(data = n_hauls,
                mapping = aes(x = SIZE_BIN, y = n_positive, linetype = TREATMENT)) +
  geom_point(data = n_hauls,
            mapping = aes(x = SIZE_BIN, y = n_positive, shape = TREATMENT)) +
  scale_color_manual(values = c("grey70", "grey20")) +
  scale_shape(solid = FALSE) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(name = "Encounters (# hauls)") +
  theme_bw() +
  theme(legend.title = element_blank())

lgcp_title <- ggdraw() + 
  draw_label(paste0(common_name, " selectivity ratio (Thygesen method)"), fontface = 'bold', x = 0.5, hjust = 0.5, size = 16)

p_lgcp <- 
  cowplot::plot_grid(
    lgcp_title,
    cowplot::plot_grid(
      p_encounters + theme(legend.position = "inside", legend.position.inside = c(0.15, 0.82),
                           legend.background = element_blank()),
      lgcp_plots$p_cpue + theme(legend.position = "none"),
      lgcp_plots$p_fit,
      nrow = 1),
    nrow = 2, rel_heights = c(0.1, 0.90)
  )

png(filename = here::here("analysis", "15_30", "plots", "selectivity_ratios", paste0(spp_code, "_lgcp_selectivity.png")),
    width = 8,
    height = 4,
    units = "in",
    res = 300)
print(p_lgcp)
dev.off()


png(filename = here::here("analysis", "15_30", "plots", "selectivity_ratios", paste0(spp_code, "_encounters_bootstrap.png")),
    width = 6,
    height = 4,
    units = "in",
    res = 300)
print(
  cowplot::plot_grid(
    draw_label(common_name, fontface = 'bold', x = 0.5, hjust = 0.5, size = 16),
    cowplot::plot_grid(
      p_encounters + theme(legend.position = "inside", legend.position.inside = c(0.15, 0.82),
                           legend.background = element_blank()),
      lgcp_plots$p_cpue + theme(legend.position = "none"),
      nrow = 1),
    nrow = 2, rel_heights = c(0.1, 0.90)
  )
)
dev.off()


# Fit to bootstrap samples
lgcp_boot <-
  lapply(
    X = boot_dat$long,
    use_logsd = use_logsd,
    FUN = function(x, use_logsd) {

      min_size <- min(x$SIZE_BIN)

      transformed <-
        x |>
        dplyr::mutate(TOTAL_COUNT = round(FREQUENCY * SAMPLING_FACTOR)) |>
        dplyr::select(-SAMPLING_FACTOR, -FREQUENCY, -ORIGNAL_MATCHUP) |>
        dplyr::group_by(HAULJOIN, MATCHUP, TREATMENT, AREA_SWEPT_KM2, SIZE_BIN) |>
        dplyr::summarise(TOTAL_COUNT = sum(TOTAL_COUNT), .groups = "keep") |>
        dplyr::ungroup() |>
        dplyr::arrange(SIZE_BIN) |>
        tidyr::pivot_wider(values_from = "TOTAL_COUNT", names_from = "SIZE_BIN", values_fill = 0);

      start_index <- which(names(transformed) == min_size)

      input <-
        list(
          N = as.matrix(transformed[, start_index:ncol(transformed)]),
          SweptArea = transformed$AREA_SWEPT_KM2,
          group = factor(transformed$MATCHUP),
          Gear = factor(transformed$TREATMENT),
          Lvec = as.numeric(names(transformed)[start_index:ncol(transformed)])
        )

      if(use_logsd) {
        fit <- gearcalib_fit(d = input, model = "poisson")
      } else {
        fit <- gearcalib_fit(d = input, model = "poisson", logsdGearRW = -10)
      }

      output <-
        data.frame(SIZE_BIN = fit$d$Lvec,
                   fit = fit$est,
                   sd_fit = fit$sd,
                   s12 = exp(fit$est))

      return(output)

    }
  )

lgcp_bootstrap_quantiles <-
  do.call(what = rbind, arg = lgcp_boot) |>
  dplyr::mutate(
    SPECIES_CODE = spp_code,
    common_name = sratio::species_code_label(x = spp_code, type = "common_name")
  ) |>
  dplyr::group_by(SIZE_BIN, SPECIES_CODE, common_name) |>
  dplyr::summarise(sratio_q025 = quantile(s12, 0.025, na.rm = TRUE),
                   sratio_q250 = quantile(s12, 0.25, na.rm = TRUE),
                   sratio_q500 = quantile(s12, 0.5, na.rm = TRUE),
                   sratio_q750 = quantile(s12, 0.75, na.rm = TRUE),
                   sratio_q975 = quantile(s12, 0.975, na.rm = TRUE)) |>
  dplyr::mutate(method = "LGCP", agg_level = "haul")


p_lgcp <-
  ggplot() +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(
    data = survey_level_bootstrap,
    mapping = aes(x = SIZE_BIN,
                  ymin = s_q025,
                  ymax = s_q975),
    width = 0
  ) +
  geom_point(data = survey_level_bootstrap,
            mapping = aes(x = SIZE_BIN, y = s_q500)
  ) +
  geom_ribbon(
    data = lgcp_bootstrap_quantiles,
    mapping = aes(x = SIZE_BIN,
                  ymin = sratio_q025,
                  ymax = sratio_q975),
    alpha = 0.2
  ) +
  geom_path(
    data = lgcp_bootstrap_quantiles,
    mapping = aes(x = SIZE_BIN,
                  y = sratio_q500)
  ) +
  geom_vline(xintercept = 10) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(
    name = "Relative selectivity",
    expand = c(0, 0),
    limits = c(0, 2),
    oob = scales::squish_infinite
  ) +
  theme_bw()

# Kotwicki selectivity ratio on haul-level data ----------------------------------------------------

gam_knots <- length(unique(sratio_dat$SIZE_BIN))-4

if(gam_knots > 10) {
  gam_knots <- 8
}

# Fewer knots for Alaska skate and red king crab
if(spp_code %in% c(471, 69322)) {
  gam_knots <- 5
}

# Binomial model matchup-level cross-validation
output_binomial <-
  sratio_cv(
    model_type = "binomial",
    count1 = sratio_dat$N_30,
    count2 = sratio_dat$N_15,
    effort1 = sratio_dat$AREA_SWEPT_KM2_30,
    effort2 = sratio_dat$AREA_SWEPT_KM2_15,
    sampling_factor1 = sratio_dat$SAMPLING_FACTOR_30,
    sampling_factor2 = sratio_dat$SAMPLING_FACTOR_15,
    size = sratio_dat$SIZE_BIN,
    block = sratio_dat$MATCHUP,
    k = gam_knots,
    n_cores = 4,
    scale_method = "sv",
    sratio_type = "absolute",
    obs_weight_control =
      list(method = "count",
           max_count = 50,
           residual_type = "absolute",
           normalize_weight = FALSE)
  )

sratio_binomial_haul <- output_binomial$cv
sratio_binomial_haul$model <- output_binomial$model_settings$model_type
sratio_binomial_haul$k <- output_binomial$model_settings$k
sratio_binomial_haul$obs_weight_method <- output_binomial$model_settings$obs_weight_control$method
sratio_binomial_haul$obs_weight_max_count <- output_binomial$model_settings$obs_weight_control$max_count
sratio_binomial_haul$obs_weight_residual_type <- output_binomial$model_settings$obs_weight_control$residual_type
sratio_binomial_haul$obs_weight_normalize_weight <- output_binomial$model_settings$obs_weight_control$normalize_weight

# Beta regression model matchup-level cross-validation
output_beta <-
  sratio_cv(
    model_type = "beta",
    count1 = sratio_dat$N_30,
    count2 = sratio_dat$N_15,
    effort1 = sratio_dat$AREA_SWEPT_KM2_30,
    effort2 = sratio_dat$AREA_SWEPT_KM2_15,
    sampling_factor1 = sratio_dat$SAMPLING_FACTOR_30,
    sampling_factor2 = sratio_dat$SAMPLING_FACTOR_15,
    size = sratio_dat$SIZE_BIN,
    block = sratio_dat$MATCHUP,
    k = gam_knots,
    n_cores = 4,
    scale_method = "sv",
    sratio_type = "absolute",
    obs_weight_control =
      list(method = "count",
           max_count = 50,
           residual_type = "none",
           normalize_weight = FALSE)
  )

sratio_beta_haul <- output_beta$cv
sratio_beta_haul$model <- output_beta$model_settings$model_type
sratio_beta_haul$k <- output_beta$model_settings$k
sratio_beta_haul$obs_weight_method <- output_beta$model_settings$obs_weight_control$method
sratio_beta_haul$obs_weight_max_count <- output_beta$model_settings$obs_weight_control$max_count
sratio_beta_haul$obs_weight_residual_type <- output_beta$model_settings$obs_weight_control$residual_type
sratio_beta_haul$obs_weight_normalize_weight <- output_beta$model_settings$obs_weight_control$normalize_weight

# Rename columns to match inputs
sratio_haul <-
  dplyr::bind_rows(sratio_binomial_haul, sratio_beta_haul) |>
  dplyr::mutate(SPECIES_CODE = spp_code) |>
  dplyr::select(
    model,
    obs_weight_method,
    obs_weight_max_count,
    obs_weight_residual_type,
    obs_weight_normalize_weight,
    k,
    SPECIES_CODE,
    SIZE_BIN = size,
    MATCHUP = block,
    N_30 = count1,
    N_15 = count2,
    SAMPLING_FACTOR_30 = sampling_factor1,
    SAMPLING_FACTOR_15 = sampling_factor2,
    AREA_SWEPT_KM2_30 = effort1,
    AREA_SWEPT_KM2_15 = effort2,
    p,
    s,
    p_fit,
    s_fit
  )

# Calculate root mean square error for proportions
sratio_haul_rmse <-
  sratio_haul |>
  dplyr::group_by(
    SPECIES_CODE,
    model,
    k,
    obs_weight_method,
    obs_weight_max_count,
    obs_weight_residual_type,
    obs_weight_normalize_weight
  ) |>
  dplyr::summarise(
    rmse = sqrt(mean((p_fit - p)^2))
  )

sratio_haul_rmse$best <- sratio_haul_rmse$rmse == min(sratio_haul_rmse$rmse)

# Fit best model to bootstrap samples
sratio_haul_bootstrap_fit <-
  sratio::sratio_fit_bootstrap(
    x = boot_dat$long,
    treatment_order = c(30, 15),
    size_col = "SIZE_BIN",
    block_col = "MATCHUP",
    treatment_col = "TREATMENT",
    count_col = "FREQUENCY",
    effort_col = "AREA_SWEPT_KM2",
    sampling_factor_col = "SAMPLING_FACTOR",
    gam_family = sratio_haul_rmse$model[sratio_haul_rmse$best],
    obs_weight_control =
      list(method = "count",
           max_count = 50,
           residual_type = "none",
           normalize_weight = FALSE),
    k = gam_knots,
    scale_method = "sv",
    sratio_type = "absolute",
    n_cores = 4
  ) |>
  dplyr::mutate(
    SPECIES_CODE = spp_code,
    common_name = sratio::species_code_label(x = SPECIES_CODE, type = "common_name")
      )

sratio_haul_bootstrap_quantiles <-
  sratio_haul_bootstrap_fit |>
  dplyr::group_by(SIZE_BIN, SPECIES_CODE, common_name) |>
  dplyr::summarise(p_q025 = quantile(p12, 0.025),
                   p_q250 = quantile(p12, 0.25),
                   p_q500 = quantile(p12, 0.5),
                   p_q750 = quantile(p12, 0.75),
                   p_q975 = quantile(p12, 0.975),
                   sratio_q025 = quantile(s12, 0.025),
                   sratio_q250 = quantile(s12, 0.25),
                   sratio_q500 = quantile(s12, 0.5),
                   sratio_q750 = quantile(s12, 0.75),
                   sratio_q975 = quantile(s12, 0.975)) |>
  dplyr::mutate(method = "Selectivity ratio", agg_level = "haul")

# Make plots of catch ratio and selectivity ratio
hist_df <-
  sratio_haul |>
  dplyr::mutate(MATCHUP = as.numeric(as.character(MATCHUP))) |>
  dplyr::inner_join(sratio::data_1530$haul |>
                      dplyr::select(MATCHUP, YEAR) |>
                      unique(),
                    by = "MATCHUP") |>
  dplyr::select(MATCHUP, SIZE_BIN, YEAR) |>
  unique()

plot_obs_histogram <-
  ggplot() +
  geom_histogram(data = hist_df,
                 mapping = aes(x = SIZE_BIN, fill = factor(YEAR)),
                 bins = length(unique(hist_df $SIZE_BIN))-1) +
  scale_x_continuous(name = xlab, expand = c(0,0)) +
  scale_y_continuous(name = "Pairs (#)") +
  scale_fill_manual(values = year_colors) +
  theme_bw() +
  theme(legend.position = c(0.17,0.87),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 6.5),
        legend.key.height = unit(2, units = "mm"),
        legend.key.width = unit(4, units = "mm"))

### ADD point colors for years!!!!!
plot_pratio <-
  ggplot() +
  geom_point(data = sratio_haul,
             mapping = aes(x = SIZE_BIN, y = p),
             size = rel(0.3),
             alpha = 0.5) +
  geom_ribbon(data = sratio_haul_bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            ymin = p_q025,
                            max = p_q975),
              alpha = 0.2,
              fill = "grey20") +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = p_q250),
            linetype = 3) +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = p_q750),
            linetype = 3) +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = p_q500)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(
    name = "Catch comparison rate",
                     limits = c(0, 1.05),
                     expand = c(0, 0),
                     oob = scales::squish_infinite
    ) +
  scale_color_manual(values = year_colors) +
  theme_bw()

### ADD point colors for years!!!!!
plot_sratio <-
  ggplot() +
  geom_point(data = sratio_haul,
             mapping = aes(x = SIZE_BIN, y = p/(1-p)),
             size = rel(0.3),
             alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(data = sratio_haul_bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            ymin = sratio_q025,
                            max = sratio_q975),
              alpha = 0.2,
              fill = "grey20") +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q250),
            linetype = 3) +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q750),
            linetype = 3) +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q500)) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(
    name = "Relative selectivity",
                expand = c(0, 0),
                limits = c(0, 2),
                oob = scales::squish_infinite
    ) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw()

plot_sratio_boot <-
  ggplot() +
  geom_errorbar(
    data = survey_level_bootstrap,
    mapping = aes(x = SIZE_BIN,
                  ymin = s_q025,
                  ymax = s_q975),
    width = 0,
    linewidth = rel(0.1)
  ) +
  geom_point(data = survey_level_bootstrap,
             mapping = aes(x = SIZE_BIN, y = s_q500),
             size = rel(0.3),
  ) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_ribbon(data = sratio_haul_bootstrap_quantiles,
              mapping = aes(x = SIZE_BIN,
                            ymin = sratio_q025,
                            max = sratio_q975),
              alpha = 0.2,
              fill = "grey20") +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q250),
            linetype = 3) +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q750),
            linetype = 3) +
  geom_path(data = sratio_haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN,
                          y = sratio_q500)) +
  scale_x_continuous(name = xlab) +
  scale_y_continuous(
    name = "Relative selectivity",
    expand = c(0, 0),
    limits = c(0, 2),
    oob = scales::squish_infinite
  ) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw()

# Write plots to file
ragg::agg_png(file = here::here("analysis", "15_30",
                                "plots", "sratio_fit", paste0(spp_code, "_sratio_three_panel_v1.png")),
              width = 169, height = 70, units = "mm", res = 300)
print(
  cowplot::plot_grid(
    plot_obs_histogram,
    plot_pratio,
    plot_sratio,
    nrow = 1,
    labels = LETTERS[1:3]
  )
)
dev.off()

ragg::agg_png(file = here::here("analysis", "15_30",
                                "plots", "sratio_fit", paste0(spp_code, "_sratio_three_panel_v2.png")),
              width = 169, height = 70, units = "mm", res = 300)
print(
  cowplot::plot_grid(
    plot_obs_histogram,
    plot_pratio,
    plot_sratio_boot,
    nrow = 1,
    labels = LETTERS[1:3]
  )
)
dev.off()

ragg::agg_png(file = here::here("analysis", "15_30",
                                "plots", "sratio_fit",
                                paste0(spp_code, "_sratio_two_panel.png")),
              width = 104, height = 70, units = "mm", res = 300)
print(
  cowplot::plot_grid(
    plot_obs_histogram,
    plot_sratio_boot,
    nrow = 1,
    labels = LETTERS[1:3]
  )
)
dev.off()

# selfisher (Brooks et al. (2022) on haul-level data -----------------------------------------------

selfisher_haul_dat <-
  sratio_dat |>
  dplyr::mutate( # offset_q based on sampling factor and area swept
    offset_q = AREA_SWEPT_KM2_30/AREA_SWEPT_KM2_15 * SAMPLING_FACTOR_15 / SAMPLING_FACTOR_30
    )

mean_size <- mean(rep(selfisher_haul_dat$SIZE_BIN, selfisher_haul_dat$N_TOTAL))
var_size <- var(rep(selfisher_haul_dat$SIZE_BIN, selfisher_haul_dat$N_TOTAL))

selfisher_haul_dat$scaled_size <- (selfisher_haul_dat$SIZE_BIN-mean_size)/sqrt(var_size)

selfisher_haul_mod <-
  selfisher::selfisher(
    p12 ~ offset(log(offset_q)) + bs(scaled_size, df = gam_knots) + (1 | MATCHUP),
    data = selfisher_haul_dat,
    total = N_TOTAL,
    haul = MATCHUP,
    psplit = FALSE
)

# Bootstrap estimate confidence intervals (REPLACE WITH PRE-DRAWN BOOTSTRAP SAMPLES)

sf_haul_boot <-
  lapply(X = boot_dat$wide,
       gam_knots = gam_knots,
       mean_size = mean_size,
       var_size = var_size,
       FUN =
         function(x, gam_knots, mean_size, var_size) {
           boot_sel <-
             x |>
             dplyr::mutate(
               scaled_size = (SIZE_BIN-mean_size)/sqrt(var_size),
               offset_q = AREA_SWEPT_KM2_30/AREA_SWEPT_KM2_15 * SAMPLING_FACTOR_15 / SAMPLING_FACTOR_30,
               N_TOTAL = N_15+N_30,
               p12 = N_30/(N_TOTAL)
             )

           fit <-
             data.frame(
               SIZE_BIN = min(boot_sel$SIZE_BIN):max(boot_sel$SIZE_BIN),
               offset_q = 1,
               N_TOTAL = 1,
               MATCHUP = 999
             ) |> # Dummy matchup
             dplyr::mutate(scaled_size = (SIZE_BIN-mean_size)/sqrt(var_size))

           sf_haul_boot_mod <-
             selfisher::selfisher(
               p12 ~ offset(log(offset_q)) + bs(scaled_size, df = gam_knots) + (1 | MATCHUP),
               data = boot_sel,
               total = N_TOTAL,
               haul = MATCHUP,
               psplit = FALSE
             )

           fit$fit <-
             predict(
               sf_haul_boot_mod,
               newdata = fit,
               type = "ratio",
               allow.new.levels = TRUE) # Fixed effects only

           return(fit)

         }

)

sf_haul_bootstrap_quantiles <-
  do.call(what = rbind, args = sf_haul_boot) |>
  dplyr::group_by(SIZE_BIN) |>
  dplyr::summarise(
    sratio_q025 = quantile(fit, 0.025),
    sratio_q250 = quantile(fit, 0.25),
    sratio_q500 = quantile(fit, 0.5),
    sratio_q750 = quantile(fit, 0.75),
    sratio_q975 = quantile(fit, 0.975)
  ) |>
  dplyr::mutate(method = "selfisher-bs", agg_level = "haul")

p_selfisher_haul <-
  ggplot()+
  geom_errorbar(
    data = survey_level_bootstrap,
    mapping = aes(x = SIZE_BIN,
                  ymin = s_q025,
                  ymax = s_q975),
    width = 0
  ) +
  geom_point(data = survey_level_bootstrap,
             mapping = aes(x = SIZE_BIN, y = s_q500)
  ) +
  geom_ribbon(
    data = sf_haul_bootstrap_quantiles,
    mapping = aes(x = SIZE_BIN, ymin = sratio_q025, ymax = sratio_q975),
    alpha = 0.2
  ) +
  geom_line(
    data = sf_haul_bootstrap_quantiles,
    mapping = aes(x = SIZE_BIN, y = sratio_q500)
  ) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(
    name = "Relative Selectivity",
    expand = c(0, 0),
    limits = c(0, 2),
    oob = scales::squish_infinite
  ) +
  scale_x_continuous(name = xlab) +
  theme_bw()

# Miller binomial on haul-level data ---------------------------------------------------------------

# INCOMPLETE

# Miller betabinomial on haul-level data -----------------------------------------------------------

# INCOMPLETE

# Plot results of all haul-level estimation methods ------------------------------------------------

haul_bootstrap_quantiles <-
  dplyr::bind_rows(
    lgcp_bootstrap_quantiles,
    sf_haul_bootstrap_quantiles,
    sratio_haul_bootstrap_quantiles
  )

p_haul_methods <-
  ggplot() +
  geom_errorbar(
    data = survey_level_bootstrap,
    mapping = aes(x = SIZE_BIN,
                  ymin = s_q025,
                  ymax = s_q975),
    width = 0,
    color = "grey30",
    alpha = 0.7,
    linewidth = rel(0.2)
  ) +
  geom_point(data = survey_level_bootstrap,
             mapping = aes(x = SIZE_BIN, y = s_q500),
             color = "grey30",
             size = rel(0.4)
  ) +
  geom_ribbon(data = haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN, ymin = sratio_q025, ymax = sratio_q975, fill = method),
            alpha = 0.3) +
  geom_path(data = haul_bootstrap_quantiles,
            mapping = aes(x = SIZE_BIN, y = sratio_q500, color = method),
            linewidth = rel(0.5)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = rel(0.4)) +
  scale_y_continuous(
    name = "Relative Selectivity",
    expand = c(0, 0),
    limits = c(0, 2),
    oob = scales::squish_infinite
  ) +
  scale_x_continuous(name = xlab) +
  scale_color_manual(name = NULL, values = c("#4C413FFF", "#278B9AFF", "#E75B64FF", "#DE7862FF", "#D8AF39FF")) +
  scale_fill_manual(name = NULL, values = c("#4C413FFF", "#278B9AFF", "#E75B64FF", "#DE7862FF", "#D8AF39FF")) +
  facet_wrap(~method) +
  theme_bw() +
  theme(legend.position = "none")

ragg::agg_png(file = here::here("analysis", "15_30",
                                "plots", "sratio_fit", paste0(spp_code, "_haul_sratio_methods.png")),
              width = 169, height = 70, units = "mm", res = 300)
print(p_haul_methods)
dev.off()

# Binomial and beta regression selectivity ratio on pooled data ------------------------------------
# K-fold cross-validation for model selection?

sratio_pooled_dat <-
  sratio_dat |>
  dplyr::group_by(SPECIES_CODE, SIZE_BIN) |>
  dplyr::summarise(
    N_OBS = as.integer(sum(N_15+N_30)),
    SUM_N_15 = sum(N_15),
    SUM_N_30 = sum(N_30),
    N_TOTAL = sum(N_15) + sum(N_30),
    SUM_FREQ_15 = sum(N_15*SAMPLING_FACTOR_15),
    SUM_FREQ_30 = sum(N_30*SAMPLING_FACTOR_30),
    SUM_AREA_SWEPT_KM2_15 = sum(AREA_SWEPT_KM2_15),
    SUM_AREA_SWEPT_KM2_30 = sum(AREA_SWEPT_KM2_30),
    offset_q = sum(AREA_SWEPT_KM2_30)/sum(AREA_SWEPT_KM2_15), # Brooks offset
    .groups = "keep"
  ) |>
  dplyr::mutate(
    TOTAL_CPUE_15 = SUM_FREQ_15/SUM_AREA_SWEPT_KM2_15,
    TOTAL_CPUE_30 = SUM_FREQ_30/SUM_AREA_SWEPT_KM2_30,
    p12_sratio = TOTAL_CPUE_30 / (TOTAL_CPUE_30 + TOTAL_CPUE_15), # Kotwicki
    s12_sratio = p12_sratio/(1-p12_sratio), # Kotwicki
    p12_sf = SUM_FREQ_30/(SUM_FREQ_30+SUM_FREQ_15), # Brooks
    s12_sf = p12_sf/(1-p12_sf) # Brooks
  )

# Data frame for generating predictions
sratio_pooled_fit <-
  data.frame(SIZE_BIN = min(sratio_pooled_dat$SIZE_BIN):max(sratio_pooled_dat$SIZE_BIN))

sratio_pooled_bin_gam  <-
  mgcv::gam(
    formula = p12_sratio ~ s(SIZE_BIN, bs = "tp", k = gam_knots),
    data = sratio_pooled_dat,
    weights = N_OBS,
    family = binomial(link = "logit")
  )

# Note: no Smithson and Verkulien (2006) transformation
sratio_pooled_beta_gam  <-
  mgcv::gam(
    formula = p12_sratio ~ s(SIZE_BIN, bs = "tp", k = gam_knots),
    data = sratio_pooled_dat,
    weights = N_OBS,
    family = betar(link = "logit")
  )

sratio_bin_pooled_fit <-
  sratio_pooled_fit |>
  cbind(
    predict(
      object = sratio_pooled_bin_gam,
      newdata = sratio_pooled_fit,
      type = "link",
      se.fit = TRUE) |>
      data.frame()
  ) |>
  dplyr::rename(logit_p12 = fit, logit_se_p12 = se.fit) |>
  dplyr::mutate(
    p12_fit = sratio::inv_logit(logit_p12),
    p12_upr = sratio::inv_logit(logit_p12 + 2 * logit_se_p12),
    p12_lwr = sratio::inv_logit(logit_p12 - 2 * logit_se_p12),
    s12_fit = p12_fit/(1 - p12_fit),
    s12_upr = p12_upr/(1 - p12_upr),
    s12_lwr = p12_lwr/(1 - p12_lwr),
    method = "SR binomial",
    agg_level = "pooled"
  )

sratio_beta_pooled_fit <-
  sratio_pooled_fit |>
  cbind(
    predict(
      object = sratio_pooled_beta_gam,
      newdata = sratio_pooled_fit,
      type = "link",
      se.fit = TRUE) |>
      data.frame()
  ) |>
  dplyr::rename(logit_p12 = fit, logit_se_p12 = se.fit) |>
  dplyr::mutate(
    p12_fit = sratio::inv_logit(logit_p12),
    p12_upr = sratio::inv_logit(logit_p12 + 2 * logit_se_p12),
    p12_lwr = sratio::inv_logit(logit_p12 - 2 * logit_se_p12),
    s12_fit = p12_fit/(1 - p12_fit),
    s12_upr = p12_upr/(1 - p12_upr),
    s12_lwr = p12_lwr/(1 - p12_lwr),
    method = "SR beta",
    agg_level = "pooled"
  )


# selfisher on pooled data -------------------------------------------------------------------------

selfisher_pooled_mod <-
  selfisher::selfisher(
    p12_sf ~ offset(log(offset_q)) + bs(SIZE_BIN, df = gam_knots),
    data = sratio_pooled_dat,
    total = N_TOTAL,
    psplit = FALSE
  )

selfisher_pooled_fit <-
  sratio_pooled_fit |>
  cbind(
    predict(
      object = selfisher_pooled_mod,
      newdata = cbind(sratio_pooled_fit, "offset_q" = 1, "N_TOTAL" = 1),
      type = "link",
      se.fit = TRUE) |>
      data.frame()
  ) |>
  dplyr::rename(logit_p12 = fit, logit_se_p12 = se.fit) |>
  dplyr::mutate(
    p12_fit = sratio::inv_logit(logit_p12),
    p12_upr = sratio::inv_logit(logit_p12 + 2 * logit_se_p12),
    p12_lwr = sratio::inv_logit(logit_p12 - 2 * logit_se_p12),
    s12_fit = p12_fit/(1 - p12_fit),
    s12_upr = p12_upr/(1 - p12_upr),
    s12_lwr = p12_lwr/(1 - p12_lwr),
    method = "selfisher",
    agg_level = "pooled"
  )

# Webster et al. (2020) on aggregate data ----------------------------------------------------------

# INCOMPLETE

# Regression w/ Tweedie on aggregate data ----------------------------------------------------------

# INCOMPLETE

pooled_fit <-
  dplyr::bind_rows(
    sratio_beta_pooled_fit,
    sratio_bin_pooled_fit,
    selfisher_pooled_fit
  )

# Plot aggregate data

p_pooled_methods <-
  ggplot() +
  geom_errorbar(
    data = survey_level_bootstrap,
    mapping = aes(x = SIZE_BIN,
                  ymin = s_q025,
                  ymax = s_q975),
    width = 0,
    color = "grey30",
    alpha = 0.7,
    linewidth = rel(0.2)
  ) +
  geom_point(data = survey_level_bootstrap,
             mapping = aes(x = SIZE_BIN, y = s_q500),
             color = "grey30",
             size = rel(0.4)
  ) +
  geom_ribbon(
    data = pooled_fit,
    mapping =
      aes(x = SIZE_BIN, ymin = s12_lwr, ymax = s12_upr, fill = method),
    alpha = 0.2
  ) +
  geom_path(
    data = pooled_fit,
    mapping = aes(x = SIZE_BIN,
                  y = s12_fit,
                  color = method)
  ) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(name = "Relative selectivity", limits = c(0, 2),
                     expand = c(0,0),
                     oob = scales::oob_squish_infinite
  ) +
  facet_wrap(~method) +
  scale_y_continuous(
    name = "Relative Selectivity",
    expand = c(0, 0),
    limits = c(0, 2),
    oob = scales::squish_infinite
  ) +
  scale_color_manual(name = NULL, values = c("#4C413FFF", "#278B9AFF", "#E75B64FF", "#DE7862FF", "#D8AF39FF")) +
  scale_fill_manual(name = NULL, values = c("#4C413FFF", "#278B9AFF", "#E75B64FF", "#DE7862FF", "#D8AF39FF")) +
  scale_x_continuous(name = xlab) +
  theme_bw() +
  theme(legend.position = "none")

ragg::agg_png(file = here::here("analysis", "15_30",
                                "plots", "sratio_fit", paste0(spp_code, "_pooled_sratio_methods.png")),
              width = 169, height = 70, units = "mm", res = 300)
print(p_pooled_methods)
dev.off()

