library(sratio)

pratio_samples <- readRDS(here::here("analysis", "15_30", "output", "pratio_samples.rds"))

pratio_samples$FIT_N_30 = sratio::sratio_predict_count(est_count = "count1",
                                                       count1 = pratio_samples$N_30,
                                                       count2 = pratio_samples$N_15,
                                                       sampling_factor1 = pratio_samples$SAMPLING_FACTOR_30,
                                                       sampling_factor2 = pratio_samples$SAMPLING_FACTOR_15,
                                                       effort1 = pratio_samples$AREA_SWEPT_KM2_30,
                                                       effort2 = pratio_samples$AREA_SWEPT_KM2_15,
                                                       s12 = pratio_samples$s_fit)

pratio_samples$FIT_N_15 = sratio::sratio_predict_count(est_count = "count2",
                                                       count1 = pratio_samples$N_30,
                                                       count2 = pratio_samples$N_15,
                                                       sampling_factor1 = pratio_samples$SAMPLING_FACTOR_30,
                                                       sampling_factor2 = pratio_samples$SAMPLING_FACTOR_15,
                                                       effort1 = pratio_samples$AREA_SWEPT_KM2_30,
                                                       effort2 = pratio_samples$AREA_SWEPT_KM2_15,
                                                       s12 = pratio_samples$s_fit)

species_code <- 21720

sel_species <- pratio_samples |>
  dplyr::filter(SPECIES_CODE == species_code) |>
  dplyr::mutate(DIFF_FIT_N_30 = N_30 - FIT_N_30)

test <- sel_species |>
  dplyr::group_by(SIZE_BIN, obs_weight_method, model) |>
  dplyr:::summarise(MEAN_DIFF_FIT_N_30 = mean(DIFF_FIT_N_30),
                    MEDIAN_DIFF_FIT_N_30 = median(DIFF_FIT_N_30))

test <- sel_species |>
  dplyr::group_by(obs_weight_method, model, SIZE_BIN) |>
  dplyr:::summarise(TOTAL_FIT_N_30 = sum(FIT_N_30),
                    TOTAL_N_30 = sum(N_30))

prediction_error <- sel_species |>
  dplyr::group_by(obs_weight_method, model, SIZE_BIN) |>
  dplyr:::summarise(TOTAL_FIT_N_30 = sum(FIT_N_30),
                    TOTAL_N_30 = sum(N_30)) |>
  dplyr::mutate(P_TOTAL_N_30 = TOTAL_FIT_N_30/TOTAL_N_30)

ggplot() + 
  geom_path(data = prediction_error,
             mapping = aes(x = SIZE_BIN, 
                           y = TOTAL_N_30-TOTAL_FIT_N_30, color = model)) +
  scale_y_continuous(name = "Observed (30) - Expected (30)") +
  facet_wrap(~obs_weight_method)

ggplot() + 
  geom_path(data = prediction_error,
            mapping = aes(x = SIZE_BIN, 
                          y = P_TOTAL_N_30, color = model)) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(name = "Expected (30)/Observed (30)") +
  facet_wrap(~obs_weight_method)

# ggplot() + 
#   geom_path(data = prediction_error,
#             mapping = aes(x = SIZE_BIN, y = MRE_N_30, color = model)) +
#   scale_y_continuous(name = expression(((N[30]-hat(N[30]))/N[30])), limits = c(-100, 100)) +
#   scale_x_continuous(name = sratio::species_code_label(x = species_code, type = "axis_label")) +
#   facet_wrap(~obs_weight_method)

write.csv(pratio_samples, here::here("pratio_samples.csv"), row.names = FALSE)
