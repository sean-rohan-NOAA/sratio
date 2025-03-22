test <- tw_results$cv_results


fit <- test |>
  dplyr::mutate(fit = fit) |>
  dplyr::select(size, fit, treatment, block) |>
  tidyr::pivot_wider(names_from = c("treatment"), values_from = "fit", names_prefix = "FIT_N_", values_fill = 0) |>
  as.data.frame()

obs_count <- test |>
  dplyr::mutate(count = count * sampling_factor) |>
  dplyr::select(size, count, treatment, block) |>
  tidyr::pivot_wider(names_from = c("treatment"), values_from = "count", names_prefix = "N_", values_fill = 0) |>
  as.data.frame()

effort <- test |>
  dplyr::select(treatment, block, effort) |>
  unique() |>
  tidyr::pivot_wider(names_from = c("treatment"), values_from = "effort", names_prefix = "EFFORT_") |>
  as.data.frame()


broh <- dplyr::full_join(fit, obs_count) |>
  dplyr::mutate(p12 = FIT_N_30/(FIT_N_15+FIT_N_30))



ggplot() +
  geom_point(data = broh,
             mapping = aes(x = size, y = p12))

ggplot() +
  geom_point(data = broh,
             mapping = aes(x = size, y = 1/(1-p12)))




test1 <- 
  MixedEffectsDemo::generate_mixed_data(
    num_groups = 10, 
    num_obs_per_group = 20, 
    intercept = 20, 
    slope = 3, 
    random_intercept_sd = 3, 
    residual_sd = 5
  )

test2 <- test1
test1$gear <- 1
test1$mult <- 1
test1$effort <- 1
test1$sampling_factor <- 1

# test2$group <- test2$group + 10
test2$gear <- 2
test2$mult <- rnorm(n = nrow(test2), mean = 2, sd = 0.25)
test2$effort <- 2
test2$sampling_factor <- 5
test2$response <- test2$response / test2$sampling_factor
test2$response <- round(test2$response)

test_comb <- rbind(test1, test2)
test_comb$response <- test_comb$response * test_comb$mult
test_comb$response <- round(test_comb$response)

test_glm <- glm(response ~ predictor + offset(log(I(1/sampling_factor))), data = test2, family = poisson(link = "log"))

exp(test_glm$coefficients)

ggplot() +
  geom_point(data = test_comb,
             mapping = aes(x = predictor, y = response, color = factor(gear)))

test_mod_dat <- dplyr::inner_join(
  test_comb |>
    dplyr::select(group, predictor, response, gear) |>
    tidyr::pivot_wider(names_from = "gear", names_prefix = "n_", values_from = "response"),
  test_comb |>
    dplyr::select(group, predictor, effort, gear) |>
    tidyr::pivot_wider(names_from = "gear", names_prefix = "effort_", values_from = "effort")
) |>
  dplyr::inner_join(
    test_comb |>
      dplyr::select(group, predictor, sampling_factor, gear) |>
      tidyr::pivot_wider(names_from = "gear", names_prefix = "sampling_factor_", values_from = "sampling_factor")
  )

test_mod_dat$group <- factor(test_mod_dat$group)

m1 <- mgcv::gam(n_2 ~ 0 + n_1 + s(group, bs = "re") + offset(I(log(effort_2/sampling_factor_2) - log(effort_1/sampling_factor_1))), 
                data = test_mod_dat, 
                family = poisson(link = "log"))

m2 <- mgcv::gam(n_2 ~ 0 + n_1 + s(group, bs = "re") + offset(I(log(effort_2/sampling_factor_2) - log(effort_1/sampling_factor_1))), 
                data = test_mod_dat, 
                family = nb(link = "log"))

m3 <- mgcv::gam(n_2 ~ 0 + n_1 + s(group, bs = "re") + offset(I(log(effort_2/sampling_factor_2) - log(effort_1/sampling_factor_1))), 
                data = test_mod_dat, 
                family = tw(link = "log"))

summary(m1)
summary(m2)
summary(m3)

exp(coef(m1))

pred_dat <- data.frame(n_1 = 25:120, 
                       effort_1 = 1, 
                       effort_2 = 2, 
                       group = 1, 
                       sampling_factor_2 = 5, 
                       sampling_factor_1 = 1)

pred_dat$poisson <- predict(m1, 
                            newdata = pred_dat, 
                            exclude = "s(group)",
                            type = "response")

pred_dat$poisson_gr <- predict(m1, 
                            newdata = pred_dat, 
                            type = "response")

pred_dat$nb <- predict(m2, 
                       newdata = pred_dat, 
                       exclude = "s(group)",
                       type = "response")

pred_dat$tw <- predict(m3, 
                       newdata = pred_dat, 
                       exclude = "s(group)",
                       type = "response")

pred_dat$poisson_gr/max(pred_dat$poisson_gr) - pred_dat$poisson/max(pred_dat$poisson)

max(pred_dat$poisson_gr)/min(pred_dat$poisson_gr)

exp(coef(m1[1]))
exp(coef(m2[1]))
exp(coef(m3[1]))

sbc_cv <- function(size,
                   block,
                   count1,
                   count2,
                   effort1,
                   effort2,
                   sampling_factor1,
                   sampling_factor2,
                   gam_family,
                   gam_formula,
                   treatment_name1, 
                   treatment_name2, 
                   k = NULL, 
                   n_cores = 1
) {
  
  sel_data <- dplyr::filter(dat_sratio, SPECIES_CODE == 21720) |>
    dplyr::mutate(MATCHUP = factor(MATCHUP)) |>
    dplyr::filter(!(N_15 == 0 & N_30 == 0))
  
  gam_formula = formula(count1 ~ 0 + count2:factor(size) + s(block, bs = "re") + offset(I(log(effort1) - log(effort2))))
  gam_family = poisson(link = "log")
  size = sel_data$SIZE_BIN
  block = sel_data$MATCHUP
  count1 = sel_data$N_30
  count2 = sel_data$N_15
  effort1 = sel_data$AREA_SWEPT_KM2_30
  effort2 = sel_data$AREA_SWEPT_KM2_15
  sampling_factor1 = sel_data$SAMPLING_FACTOR_30
  sampling_factor2 = sel_data$SAMPLING_FACTOR_15
  
  stopifnot("sbc_cv: count1 and count2 must be the same length." = length(count1) == length(count2))
  
  if(length(sampling_factor1) == 1) {
    sampling_factor1 <- rep(sampling_factor1, length(size1))
  }
  
  if(length(sampling_factor2) == 1) {
    sampling_factor2 <- rep(sampling_factor2, length(size2))
  }
  
  if(length(effort1) == 1) {
    effort1 <- rep(effort1, length(size1))
  }
  
  if(length(effort2) == 1) {
    effort2 <- rep(effort2, length(size2))
  }
  
  if(!is(block, "factor")) {
    block <- as.factor(block)
  }
  
  sbc_fit_gamm <- function(data = data,
                           gam_family = poisson(link = "log"),
                           gam_formula) {
    
    mod <- mgcv::gam(formula = gam_formula,
                     family = gam_family,
                     data = data)
    
    return(mod)
    
  }
  
  model_df <- data.frame(
    size = size,
    count1 = count1 * sampling_factor1,
    count2 = count2 * sampling_factor2,
    block = block,
    effort1 = effort1,
    effort2 = effort2
  )
  
  ggplot(data = model_df,
         mapping = aes(x = count1, y = count2)) +
    geom_point() +
    geom_smooth() +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~size)
  
  index <- as.numeric(model_df$block)
  
  # test <- sbc_fit_gamm(
  #   data = model_df,
  #   gam_family = gam_family,
  #   gam_formula = gam_formula
  # )
  
  # Setup four clusters and folds for each match-ups
  doParallel::registerDoParallel(parallel::makeCluster(n_cores))
  
  cv_results <- foreach::foreach(fold = 1:max(index), .packages = "mgcv") %dopar% {
    
    training_df <- model_df[-which(index == fold), ]
    validation_df <- model_df[which((index == fold)), ]
    
    gam_formula = formula(count1 ~ count2:factor(size))
    
    mod <- sbc_fit_gamm(data = training_df,
                        gam_formula = gam_formula,
                        gam_family = gam_family)
    
    validation_df$term <- as.numeric(
      predict.gam(object = mod, 
              newdata = validation_df, 
              type = "response",
              exclude = "s(block)")
    )
    
    training_df$term <- predict(mod, type = "terms", exclude = "s(block)")
    
    training_df$response <- predict(mod, type = "response", exclude = "s(block)")
    
    dplyr::group_by(training_df, size) |>
      dplyr::summarise(diff = sum(count1 - fit))
    
    validation_df$fit <- validation_df$count2 * 
      exp(validation_df$term) * 
      (validation_df$effort1 / validation_df$effort2) * 
      (validation_df$sampling_factor2 / validation_df$sampling_factor1)
    
    return(validation_df)
  }
  
  output_cv <- do.call("rbind", cv_results)
  
  dplyr::group_by(output_cv,
                  size) |>
    dplyr::summarise(fit_diff = sum(count1-fit))
  
  doParallel::stopImplicitCluster()
  
  return(
    list(
      cv_results = output_cv,
      model_settings = 
        list(
          gam_family = gam_family,
          k = k
        ),
      data = model_df
    )
  )
  
  
}



broh2 <- sbsc_fit_gamm(data = sel_data,
                      gam_family = poisson(link = "log"),
                      gam_formula = formula(N_30 ~ 0 + s(N_15, by = factor(SIZE_BIN)) + 
                                              s(MATCHUP, bs = "re") + 
                                              offset(I(log(AREA_SWEPT_KM2_15/SAMPLING_FACTOR_15) - log(AREA_SWEPT_KM2_30/SAMPLING_FACTOR_30)))
                                            )
                      )

test_45 <- dplyr::filter(model_df, size == 45)
test_45$count2 * test_45$effort1 / test_45$effort2 * test_45$sampling_factor2 /test_45$sampling_factor1

predict(m1, 
        newdata = data.frame(n_1 = 25:120, effort_1 = 1, effort_2 = 2, group = 1, sampling_factor_2 = 2, sampling_factor_1 = 1), 
        exclude = "s(group)",
        type = "response")
