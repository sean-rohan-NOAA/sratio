library(sratio)

size = 10:55

settings = list(
  n_sims = 50,
  n_pairs = 40, 
  selectivity_opts1 = list(type = "asymptotic",
                           begin_top = 25,
                           ln_sd1 = 10) ,
  selectivity_opts2 = list(type = "asymptotic",
                           begin_top = 35,
                           ln_sd1 = 10) ,
  gear_q1 = 1,
  gear_q2 = 0.8,
  effort1 = 1,
  effort2 = 1,
  size = size,
  abundance = round(dnorm(size, mean = 35, sd = 10) * 1e5 * (1-rnorm(length(size), mean = 0, sd = 0.1))),
  demographic_comp_pars = list(distribution = "none"),
  availability = list(distribution = "normal", mean = 0.0004, sd = 0.00008),
  size_bin_width = 4,
  n_boot_samples = 100,
  gam_knots = 5,
  min_sample_size = 10,
  scale_method = "sv",
  n_cores = 5
)

all_avail <- run_selectivity_simulation(settings = settings)

plot_all_avail <- plot_simulation_results(all_avail)


cowplot::plot_grid(
  plot_all_avail$plot_rel_bias,
  plot_all_avail$plot_within_ci,
  plot_all_avail$plot_rel_precision,
  nrow = 3, align = "v"
)
