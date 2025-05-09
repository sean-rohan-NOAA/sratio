library(TMB)
library(ggplot2)
library(sratio)
library(cowplot)

TMB::compile("gearcalib_switch.cpp")
dyn.load(dynlib("gearcalib_switch"))

source("gearcalib_nb.R")

# Load output
dat <- readRDS(file = "C:/Users/sean.rohan/Work/afsc/sratio/analysis/15_30/output/catch_at_length_1530.rds")

sel_spp <- 10210
xlab <- sratio::species_code_label(x = sel_spp)

sel_dat <-
  dat |>
  dplyr::filter(SPECIES_CODE == sel_spp) |>
  dplyr::mutate(TOTAL_COUNT = round(FREQUENCY * SAMPLING_FACTOR)) |>
  dplyr::select(-SAMPLING_FACTOR, -FREQUENCY) |>
  dplyr::arrange(SIZE_BIN) |>  
  tidyr::pivot_wider(values_from = "TOTAL_COUNT", names_from = "SIZE_BIN", values_fill = 0)

d_input <- 
  list(
    N = as.matrix(sel_dat[, 6:ncol(sel_dat)]),
    SweptArea = sel_dat$AREA_SWEPT_KM2,
    group = factor(sel_dat$MATCHUP),
    Gear = factor(sel_dat$TREATMENT),
    Lvec = as.numeric(names(sel_dat)[6:ncol(sel_dat)])
  )

# Fit Poisson model ---
fit <- gearcalib_fit(d = d_input, model = "poisson")

# Fix random walk logsd at a tiny value when the Hessian is not positive definite
if(!fit$rep$pdHess) {
  fit <- gearcalib_fit(d = d_input, model = "poisson", logsdGearRW = -10)
}

fit$rep

boot_results <- gearcalib_boot(d_input, quantiles = c(0.025,0.5,0.975), nboot = 1000)

fit_plots <- gearcalib_plot(fit = fit, boot = boot_results, add_bootquantiles = TRUE, xlab = xlab)

n_hauls <-
  dat |>
  dplyr::filter(SPECIES_CODE == sel_spp) |>
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

cowplot::plot_grid(
  p_encounters + theme(legend.position = "inside", legend.position.inside = c(0.15, 0.82),
                       legend.background = element_blank()),
  fit_plots$p_cpue + theme(legend.position = "none"),
  fit_plots$p_fit,
  nrow = 1
)

