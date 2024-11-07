# Select best selectivity ratio model

library(sratio)

dat_sratio <- readRDS(file = here::here("analysis", "15_30", "output", "n_by_treatment_1530.rds"))

sp_code <- unique(dat_sratio$SPECIES_CODE)

#temp species drop
# if(!any(use_cruises %in% c(199501, 199801))) {
#   sp_code <- sp_code[-which(sp_code == 68580)]
# }

unique_matchups <- unique(dat_sratio$MATCHUP)

# n_cores <- 4

rmse_df <- data.frame()

pratio_samples <- data.frame()

for(ii in 1:length(sp_code)) {
  
  pratio_df <- data.frame()
  
  spp_lengths <- dplyr::filter(dat_sratio, SPECIES_CODE == sp_code[ii])
  
  # Set knots based on number of length bins, but only use 5 knots for red king crab
  gam_knots <- length(unique(spp_lengths$SIZE_BIN))-4
  
  if(gam_knots > 10) {
    gam_knots <- 8
  }
  
  if(sp_code[ii] %in% c(471, 69322)) {
    gam_knots <- 5
  }
  
  # Run match-up level cross validation
  pratio_df <- sratio_cv(count1 = spp_lengths$N_30,
            count2 = spp_lengths$N_15,
            effort1 = spp_lengths$AREA_SWEPT_KM2_30,
            effort2 = spp_lengths$AREA_SWEPT_KM2_15,
            size = spp_lengths$SIZE_BIN,
            block = spp_lengths$MATCHUP,
            k = gam_knots,
            n_cores = 4,
            scale_method = "sv")
  
  pratio_df$SPECIES_CODE <- sp_code[ii]
  
  # Rename columns to match inputs
  pratio_df <- pratio_df |> 
    dplyr::rename(MATCHUP = block,
                  SIZE_BIN = size,
                  N_30 = count1,
                  N_15 = count2,
                  AREA_SWEPT_KM2_30 = effort1,
                  AREA_SWEPT_KM2_15 = effort2)
  
  pratio_samples <- dplyr::bind_rows(pratio_samples, pratio_df)

  rmse_df <- rmse_df |>
    dplyr::bind_rows(
      data.frame(model = c("binomial", "beta"), 
                 rmse = c(sqrt(mean((pratio_df$cv_fit_logit-pratio_df$p)^2)),
                          sqrt(mean((pratio_df$cv_fit_beta-pratio_df$p)^2))),
                 SPECIES_CODE = sp_code[ii],
                 gam_knots = gam_knots) |>
        dplyr::mutate(best = rmse == min(rmse))
    )
}

saveRDS(pratio_samples, file = here::here("analysis", "15_30", "output", "pratio_samples.rds"))

write.csv(rmse_df, here::here("analysis", "15_30", "output", "sratio_model_rmse.csv"), row.names = FALSE)
