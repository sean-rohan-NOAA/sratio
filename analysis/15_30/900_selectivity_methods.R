library(sratio)
library(selfisher)
library(splines)
library(cowplot)

spp_code <- 21720
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
  catch_at_size_dat|>
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
  dplyr::mutate(type = "Bootstrap")

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
              alpha = 0.5,
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
    name = expression(italic(p['L,30,15'])), 
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
              alpha = 0.5,
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
    name = expression(italic(S['L,30,15'])~(SR)), 
                expand = c(0, 0),
                limits = c(0, 2),
                oob = scales::squish_infinite
    ) +
  scale_color_tableau() +
  scale_fill_tableau() +
  theme_bw()

# Write plots to file

ragg::agg_png(file = here::here("analysis", "15_30", 
                                "plots", "sratio_fit", paste0(spp_code, "_sratio_three_panel.png")), 
              width = 169, height = 70, units = "mm", res = 300)
print(cowplot::plot_grid(plot_obs_histogram,
                         plot_pratio,
                         plot_sratio,
                         nrow = 1,
                         labels = LETTERS[1:3]))
dev.off()

ragg::agg_png(file = here::here("analysis", "15_30",  
                                "plots", "sratio_fit",
                                paste0(spp_code, "_sratio_two_panel.png")), 
              width = 104, height = 70, units = "mm", res = 300)
print(cowplot::plot_grid(plot_obs_histogram,
                         plot_sratio,
                         nrow = 1,
                         labels = LETTERS[1:3]))
dev.off()

# selfisher (Brooks et al. (2022) on haul-level data -----------------------------------------------

# Offset using the sampling factor and area swept

selfisher_haul_dat <- 
  sratio_dat |>
  dplyr::mutate(offset_q = AREA_SWEPT_KM2_30/AREA_SWEPT_KM2_15 * SAMPLING_FACTOR_15 / SAMPLING_FACTOR_30)

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

boot_fit <- vector(mode = "list", length = n_boot)

for(ii in 1:n_boot) {
  
  boot_sel <- 
    boot_dat$wide[[ii]] |>
    dplyr::mutate(
      scaled_size = (SIZE_BIN-mean_size)/sqrt(var_size),
      offset_q = AREA_SWEPT_KM2_30/AREA_SWEPT_KM2_15 * SAMPLING_FACTOR_15 / SAMPLING_FACTOR_30,
      N_TOTAL = N_15+N_30,
      p12 = N_30/(N_TOTAL)
    )
  
  boot_fit[[ii]] <- 
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
  
  boot_fit[[ii]]$fit <-
    predict(
      sf_haul_boot_mod, 
      newdata = boot_fit[[ii]], 
      type = "ratio",
      allow.new.levels = TRUE)
  
  print(ii)
  
}

selfisher_haul_fit <- 
  do.call(what = rbind, args = boot_fit) |>
  dplyr::group_by(SIZE_BIN) |>
  dplyr::summarise(
    sratio_q025 = quantile(fit, 0.025),
    sratio_q250 = quantile(fit, 0.25),
    sratio_q500 = quantile(fit, 0.5),
    sratio_q750 = quantile(fit, 0.75),
    sratio_q975 = quantile(fit, 0.975)
  )

p_selfisher_haul <-
  ggplot()+
  geom_point()+
  geom_ribbon(
    data = selfisher_haul_fit, 
    mapping = aes(x = SIZE_BIN, ymin = sratio_q025, ymax = sratio_q975), 
    alpha = 0.2
  ) +
  geom_line(
    data = selfisher_haul_fit, 
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


# Miller betabinomial on haul-level data -----------------------------------------------------------


# Binomial selectivity ratio on pooled data --------------------------------------------------------
# K-fold cross-validation for model selection?

unique_matchups <- unique(sratio_dat$MATCHUP)

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

sratio_pooled_fit <- 
  data.frame(SIZE_BIN = min(sratio_pooled_dat$SIZE_BIN):max(sratio_pooled_dat$SIZE_BIN))

sratio_pooled_fit[, c("binom_logit_p12", "binom_logit_se_p12")] <- 
  predict(object = sratio_pooled_bin_gam,
          newdata = sratio_pooled_fit,
          type = "link",
          se.fit = TRUE) |>
  data.frame()

sratio_pooled_fit[, c("beta_logit_p12", "beta_logit_se_p12")] <- 
  predict(object = sratio_pooled_beta_gam,
          newdata = sratio_pooled_fit,
          type = "link",
          se.fit = TRUE) |>
  data.frame()

sratio_pooled_fit <-
  sratio_pooled_fit |>
  dplyr::mutate(
    binom_p12_fit = sratio::inv_logit(binom_logit_p12),
    binom_p12_upr = sratio::inv_logit(binom_logit_p12 + 2 * binom_logit_se_p12),
    binom_p12_lwr = sratio::inv_logit(binom_logit_p12 - 2 * binom_logit_se_p12),
    binom_s12_fit = binom_p12_fit/(1 - binom_p12_fit),
    binom_s12_upr = binom_p12_upr/(1 - binom_p12_upr),
    binom_s12_lwr = binom_p12_lwr/(1 - binom_p12_lwr),
    beta_p12_fit = sratio::inv_logit(beta_logit_p12),
    beta_p12_upr = sratio::inv_logit(beta_logit_p12 + 2 * beta_logit_se_p12),
    beta_p12_lwr = sratio::inv_logit(beta_logit_p12 - 2 * beta_logit_se_p12),
    beta_s12_fit = beta_p12_fit/(1 - beta_p12_fit),
    beta_s12_upr = beta_p12_upr/(1 - beta_p12_upr),
    beta_s12_lwr = beta_p12_lwr/(1 - beta_p12_lwr)
  )


# selfisher on aggregate data ----------------------------------------------------------------------

selfisher_pooled_mod <- 
  selfisher::selfisher(
    p12_sf ~ offset(log(offset_q)) + bs(SIZE_BIN, df = gam_knots), 
    data = sratio_pooled_dat, 
    total = N_TOTAL,
    psplit = FALSE
  )

sratio_pooled_fit[, c("sf_logit_p12", "sf_logit_se_p12")] <- 
  predict(
    object = selfisher_pooled_mod, 
    newdata = cbind(sratio_pooled_fit, "offset_q" = 1, "N_TOTAL" = 1), 
    type = "link", 
    se.fit = TRUE
  ) |>
  as.data.frame()

predict(
  object = selfisher_pooled_mod, 
  newdata = cbind(sratio_pooled_fit, "offset_q" = 1, "N_TOTAL" = 1), 
  type = "ratio", 
  se.fit = TRUE
)

sratio_pooled_fit <-
  sratio_pooled_fit |>
  dplyr::mutate(
    sf_p12_fit = sratio::inv_logit(sf_logit_p12),
    sf_p12_upr = sratio::inv_logit(sf_logit_p12 + 2 * sf_logit_se_p12),
    sf_p12_lwr = sratio::inv_logit(sf_logit_p12 - 2 * sf_logit_se_p12),
    sf_s12_fit = sf_p12_fit/(1 - sf_p12_fit),
    sf_s12_upr = sf_p12_upr/(1 - sf_p12_upr),
    sf_s12_lwr = sf_p12_lwr/(1 - sf_p12_lwr)
  )

# Regression w/ Tweedie on aggregate data ----------------------------------------------------------



# Plot aggregate data  
ggplot() +
  geom_ribbon(
    data = sratio_pooled_fit, 
    mapping = 
      aes(x = SIZE_BIN, ymin = binom_s12_lwr, ymax = binom_s12_upr, fill = "Binomial GAM"), 
    alpha = 0.2
  ) +
  geom_path(data = sratio_pooled_fit,
            mapping = aes(x = SIZE_BIN,
                          y = binom_s12_fit,
                          color = "Binomial GAM")) +
  geom_ribbon(
    data = sratio_pooled_fit, 
    mapping = 
      aes(x = SIZE_BIN, ymin = beta_s12_lwr, ymax = beta_s12_upr, fill = "Beta GAM"), 
    alpha = 0.2
  ) +
  geom_path(data = sratio_pooled_fit,
            mapping = aes(x = SIZE_BIN,
                          y = beta_s12_fit,
                          color = "Beta GAM")) +
  geom_ribbon(
    data = sratio_pooled_fit,
    mapping =
      aes(x = SIZE_BIN, ymin = sf_s12_lwr, ymax = sf_s12_upr, fill = "selfisher spline"),
    alpha = 0.2
  ) +
  geom_path(data = sratio_pooled_fit,
            mapping = aes(x = SIZE_BIN,
                          y = sf_s12_fit,
                          color = "selfisher spline")) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(name = "Relative selectivity", limits = c(0, 2),
                     expand = c(0,0),
                     oob = scales::oob_squish_infinite
  ) +
  scale_y_continuous(
    name = "Relative Selectivity",
    expand = c(0, 0),
    limits = c(0, 2),
    oob = scales::squish_infinite
  ) +
  scale_color_manual(name = NULL, values = c("#4C413FFF", "#278B9AFF", "#E75B64FF", "#DE7862FF", "#D8AF39FF")) +
  scale_fill_manual(name = NULL, values = c("#4C413FFF", "#278B9AFF", "#E75B64FF", "#DE7862FF", "#D8AF39FF")) +
  scale_x_continuous(name = xlab) +
  theme_bw()

# Webster et al. (2020) on aggregate data ----------------------------------------------------------

