# Draw bootstrap samples for selectivity ratio and total catch models for 15/30

library(sratio)

# Treatment levels
treatments <- factor(c(15,30))

seed <- 19710
n_boot_draws <- 1000

dat <- readRDS(here::here("output", "catch_at_length_1530.rds"))

sp_code <- unique(dat$SPECIES_CODE)

for(jj in 1:length(sp_code)) {
  
  dir.create(here::here("output", sp_code[jj]))
  
  sel_dat <- dplyr::filter(dat, SPECIES_CODE == sp_code[jj]) |>
    dplyr::group_by(HAULJOIN, SPECIES_CODE, MATCHUP, SIZE_BIN, AREA_SWEPT_KM2, TREATMENT) |>
    dplyr::summarize(FREQ_EXPANDED = sum(FREQ_EXPANDED), .groups = "keep") |>
    as.data.frame()
  
  # set.seed(seed)
  # bootstrap_1 <- sratio::two_stage_bootstrap(x = sel_dat$SIZE_BIN[sel_dat$TREATMENT == treatments[1]],
  #                                            group = sel_dat$MATCHUP[sel_dat$TREATMENT == treatments[1]], 
  #                                            frequency = sel_dat$FREQ_EXPANDED[sel_dat$TREATMENT == treatments[1]],
  #                                            draws = n_boot_draws)
  # 
  # set.seed(seed)
  # bootstrap_2 <- sratio::nested_bootstrap(x = sel_dat$SIZE_BIN[sel_dat$TREATMENT == treatments[2]], 
  #                                         frequency = sel_dat$FREQ_EXPANDED[sel_dat$TREATMENT == treatments[2]], 
  #                                         group = sel_dat$MATCHUP[sel_dat$TREATMENT == treatments[2]], 
  #                                         draw_group = bootstrap_1$draws$group, 
  #                                         draw_index = bootstrap_1$draws$draw)
  # 
  # boot_samples <- vector(mode = "list", length = n_boot_draws)
  
  # for(kk in 1:n_boot_draws) {
  #   
  #   boot_samples[[kk]] <- rbind(
  #     dplyr::mutate(bootstrap_1$samples[[kk]],
  #                   TREATMENT = treatments[1]) |>
  #       dplyr::rename(SIZE_BIN = x, MATCHUP = group),
  #     dplyr::mutate(bootstrap_2$samples[[kk]],
  #                   TREATMENT = treatments[2]) |>
  #       dplyr::rename(SIZE_BIN = x, MATCHUP = group)
  #   ) |>
  #     dplyr::group_by(MATCHUP, TREATMENT, SIZE_BIN) |>
  #     dplyr::summarize(FREQ_EXPANDED = n(), .groups = "keep") |>
  #     dplyr::inner_join(dplyr::select(sel_dat,
  #                                     MATCHUP, TREATMENT, AREA_SWEPT_KM2, HAULJOIN) |>
  #                         unique(),
  #                       by = c('MATCHUP', 'TREATMENT')) |>
  #     as.data.frame()
  #   
  # }
  
  boot_samples <- sratio::two_stage_bootstrap(count1 = sel_dat$FREQ_EXPANDED[sel_dat$TREATMENT == treatments[1]],
                                              count2 = sel_dat$FREQ_EXPANDED[sel_dat$TREATMENT == treatments[2]],
                                              size1 = sel_dat$SIZE_BIN[sel_dat$TREATMENT == treatments[1]],
                                              size2 = sel_dat$SIZE_BIN[sel_dat$TREATMENT == treatments[2]],
                                              block1 = sel_dat$MATCHUP[sel_dat$TREATMENT == treatments[1]],
                                              block2 = sel_dat$MATCHUP[sel_dat$TREATMENT == treatments[2]],
                                              n_draws = 1000,
                                              seed = seed,
                                              treatment_name1 = treatments[1],
                                              treatment_name2 = treatments[2])
  
  # Rename output columns and join with effort
  for(ii in 1:length(boot_samples)) {
    names(boot_samples[[ii]]) <- c("MATCHUP", "TREATMENT", "SIZE_BIN", "FREQ_EXPANDED")
    
    boot_samples[[ii]] <- boot_samples[[ii]] |>
      dplyr::inner_join(dplyr::select(sel_dat, MATCHUP, TREATMENT, AREA_SWEPT_KM2, HAULJOIN) |>
                          unique(),
                        by = c('MATCHUP', 'TREATMENT')
      ) |>
      as.data.frame()

  }
  
  saveRDS(boot_samples, file = here::here("output", sp_code[jj], paste0("bootstrap_samples_", sp_code[jj], ".rds")))
  
}
