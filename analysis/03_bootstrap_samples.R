# Draw bootstrap samples for selectivity ratio and total catch models for 15/30

library(sratio)

# Treatment levels
treatments <- factor(c(15,30))

seed <- 19710
n_boot_draws <- 1000

dat <- readRDS(here::here("output", "catch_at_length_1530.rds"))

species_codes <- unique(dat$SPECIES_CODE)

for(jj in 1:length(species_codes)) {
  
  sel_dat <- dplyr::filter(dat, SPECIES_CODE == species_codes[jj]) |>
    dplyr::group_by(HAULJOIN, SPECIES_CODE, MATCHUP, SIZE_BIN, AREA_SWEPT_KM2, TREATMENT) |>
    dplyr::summarize(FREQ_EXPANDED = sum(FREQ_EXPANDED), .groups = "keep") |>
    as.data.frame()
  
  set.seed(seed)
  bootstrap_1 <- sratio::two_stage_bootstrap(x = sel_dat$SIZE_BIN[sel_dat$TREATMENT == treatments[1]],
                                             group = sel_dat$MATCHUP[sel_dat$TREATMENT == treatments[1]], 
                                             frequency = sel_dat$FREQ_EXPANDED[sel_dat$TREATMENT == treatments[1]],
                                             draws = n_boot_draws)
  
  set.seed(seed)
  bootstrap_2 <- sratio::nested_bootstrap(x = sel_dat$SIZE_BIN[sel_dat$TREATMENT == treatments[2]], 
                                          frequency = sel_dat$FREQ_EXPANDED[sel_dat$TREATMENT == treatments[2]], 
                                          group = sel_dat$MATCHUP[sel_dat$TREATMENT == treatments[2]], 
                                          draw_group = bootstrap_1$draws$group, 
                                          draw_index = bootstrap_1$draws$draw)
  
  boot_samples <- vector(mode = "list", length = n_boot_draws)
  
  for(kk in 1:n_boot_draws) {
    
    boot_samples[[kk]] <- rbind(
      dplyr::mutate(bootstrap_1$samples[[kk]],
                    TREATMENT = treatments[1]) |>
        dplyr::rename(SIZE_BIN = x, MATCHUP = group),
      dplyr::mutate(bootstrap_2$samples[[kk]],
                    TREATMENT = treatments[2]) |>
        dplyr::rename(SIZE_BIN = x, MATCHUP = group)
    ) |>
      dplyr::group_by(MATCHUP, TREATMENT, SIZE_BIN) |>
      dplyr::summarize(FREQ_EXPANDED = n(), .groups = "keep") |>
      dplyr::inner_join(dplyr::select(sel_dat,
                                      MATCHUP, TREATMENT, AREA_SWEPT_KM2, HAULJOIN) |>
                          unique(), by = c('MATCHUP', 'TREATMENT')) |>
      dplyr::mutate(MATCHUP = factor(MATCHUP)) |>
      as.data.frame()
    
  }
  
  saveRDS(boot_samples, file = here::here("output", species_codes[jj], paste0("bootstrap_samples_", species_codes[jj], ".rds")))
  
}