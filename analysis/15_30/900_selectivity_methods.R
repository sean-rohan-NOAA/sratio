library(ggplot2)
library(sratio)
library(cowplot)

spp_code <- 685801
xlab <- unique(sratio::species_code_label(x = spp_code))
common_name <- unique(sratio::species_code_label(x = spp_code, type = "common_name"))

# Setup directories if they don't exist ------------------------------------------------------------

dir.create(here::here("analysis", "15_30", "plots", "selectivity_ratios"), 
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(here::here("analysis", "15_30", "output", "selectivity_ratios"), 
           recursive = TRUE,
           showWarnings = FALSE)

# Data ---- ----------------------------------------------------------------------------------------
catch_at_size_dat <- readRDS(file = here::here("analysis", "15_30", "output", "catch_at_length_1530.rds")) |>
  dplyr::filter(SPECIES_CODE %in% spp_code)

sratio_dat <- 
  readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds")) |>
  dplyr::filter(SPECIES_CODE %in% spp_code)

boot_dat <- 
  readRDS(here::here("analysis", "15_30", "output", spp_code, paste0("bootstrap_samples_", spp_code, ".rds")))


# Thygesen et al. (2019) log-Gaussian Cox process model --------------------------------------------
lgcp_dat <-
  catch_at_size_dat |>
  dplyr::mutate(TOTAL_COUNT = round(FREQUENCY * SAMPLING_FACTOR)) |>
  dplyr::select(-SAMPLING_FACTOR, -FREQUENCY) |>
  dplyr::group_by(HAULJOIN, MATCHUP, TREATMENT, SPECIES_CODE, AREA_SWEPT_KM2, SIZE_BIN) |>
  dplyr::summarise(TOTAL_COUNT = sum(TOTAL_COUNT)) |>
  dplyr::arrange(SIZE_BIN) |>  
  tidyr::pivot_wider(values_from = "TOTAL_COUNT", names_from = "SIZE_BIN", values_fill = 0)

d_input <- 
  list(
    N = as.matrix(lgcp_dat[, 6:ncol(lgcp_dat)]),
    SweptArea = lgcp_dat$AREA_SWEPT_KM2,
    group = factor(lgcp_dat$MATCHUP),
    Gear = factor(lgcp_dat$TREATMENT),
    Lvec = as.numeric(names(lgcp_dat)[6:ncol(lgcp_dat)])
  )

fit <- gearcalib_fit(d = d_input, model = "poisson")

# Fix random walk logsd at a tiny value when the Hessian is not positive definite
if(!fit$rep$pdHess) {
  fit <- gearcalib_fit(d = d_input, model = "poisson", logsdGearRW = -10)
}

boot_results <- gearcalib_boot(d_input, quantiles = c(0.025,0.5,0.975), nboot = 1000)

fit_plots <- gearcalib_plot(fit = fit, boot = boot_results, add_bootquantiles = TRUE, xlab = xlab)

n_hauls <-
  dat |>
  dplyr::filter(SPECIES_CODE %in% spp_code) |>
  dplyr::group_by(TREATMENT, SIZE_BIN) |>
  dplyr::summarise(n_positive = n())

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
      fit_plots$p_cpue + theme(legend.position = "none"),
      fit_plots$p_fit,
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
sratio_bootstrap_fit <- 
  sratio::sratio_fit_bootstrap(
    x = boot_dat,
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
  )

  

# Kotwicki selectivity ratio on aggregate data -----------------------------------------------------


# selfisher on haul-level data ---------------------------------------------------------------------


# selfisher on aggregate data ----------------------------------------------------------------------


# Webster et al. (2020) ----------------------------------------------------------------------------


# Miller binomial on haul-level data ---------------------------------------------------------------


# Miller betabinomial on haul-level data -----------------------------------------------------------


# Regression w/ Tweedie on aggregate data ----------------------------------------------------------

