# Plot effort reduction scenarios
library(sratio)
library(akgfmaps)
library(ggridges)
library(cowplot)


seed <- NULL
ref_distance_fished_m <- 2800
distance_between_tows_m <- 600
# fish_density_n_km2 = c(30, 50, 100, 150, 250, 500, 1000, 2000)
grid_dim_m = c(3000, 3000)
open_boundary = TRUE
# draws = 100

results <- 
  expand.grid(
    sample = 1:100,
    fish_density_n_km2 = c(20, 50, 100, 500, 1000),
    cluster_radius_m = c(30, 400, 1000, 2000),
    cluster_density_n_km2 = c(0.5, 1),
    n_hauls = c(20, 40, 80),
    retention = 1,
    poisson_ran = NA,
    poisson_clu = NA,
    binom_ccr_ran = NA,
    binom_ccr_clu = NA,
    agg_ccr_ran = NA,
    agg_ccr_clu = NA
  )

# test version
# results <- 
#   expand.grid(
#     sample = 1:50,
#     fish_density_n_km2 = 50,
#     cluster_radius_m = 400,
#     cluster_density_n_km2 = 1,
#     n_hauls = c(20, 40), 
#     retention = 1,
#     poisson_ran = NA,
#     poisson_clu = NA,
#     binom_ccr_ran = NA,
#     binom_ccr_clu = NA,
#     agg_ccr_ran = NA,
#     agg_ccr_clu = NA
#   )


start_time <- Sys.time()

for(jj in 1:nrow(results)) {
  
  if(jj%%50 == 0) {
    print(paste0(jj, " - ", Sys.time())) 
  }
  
  # Simulate spatial distributions using Poisson and Matern clusters
  poisson_sample <- 
    try(
      sim_matern_clusters(
        fish_density_n_km2 = results$fish_density_n_km2[jj],
        grid_dim_m = grid_dim_m,
        cluster_density_n_km2 = results$cluster_density_n_km2[jj],
        cluster_radius_m = results$cluster_radius_m[jj],
        open_boundary = open_boundary,
        seed = seed, 
        draws = results$n_hauls[jj]
      ),
      silent = TRUE)
  
  
  encounters <-
    try(
      data.frame(
        draw = 1:results$n_hauls[jj],
        enc_ran_30 = 
          unlist(
            lapply(
              X = poisson_sample[['random_points']], 
              FUN = sim_encounters, 
              origin_m = c(0, grid_dim_m[2]/2), 
              effort_m = c(ref_distance_fished_m, 17)
            )
          ),
        enc_ran_15 = 
          unlist(
            lapply(
              X = poisson_sample[['random_points']], 
              FUN = sim_encounters, 
              origin_m = c(0, grid_dim_m[2]/2 + distance_between_tows_m), 
              effort_m = c(ref_distance_fished_m/2, 17)
            )
          ),
        enc_clu_30 = 
          unlist(
            lapply(
              X = poisson_sample[['cluster_points']], 
              FUN = sim_encounters, 
              origin_m = c(0, grid_dim_m[2]/2), 
              effort_m = c(ref_distance_fished_m, 17)
            )
          ),
        enc_clu_15 = 
          unlist(
            lapply(
              X = poisson_sample[['cluster_points']], 
              FUN = sim_encounters, 
              origin_m = c(0, grid_dim_m[2]/2 + distance_between_tows_m), 
              effort_m = c(ref_distance_fished_m/2, 17)
            )
          )
      ),
      silent = TRUE)
  
  if(is(encounters, "try-error")) {
    next
  }
  
  # Simulate retention ----
  set.seed(seed)
  
  for(ii in 2:ncol(encounters)) {
    
    encounters <- cbind(
      encounters, 
      mapply(function(n) rbinom(1, n, results$retention[jj]), encounters[, ii])
    )
    
    names(encounters)[ncol(encounters)] <- gsub(x = names(encounters)[ii], pattern = "enc", replacement = "ret")
    
  }
  
  
  encounters$total_ret_ran <- encounters$ret_ran_15 + encounters$ret_ran_30
  encounters$total_ret_clu <- encounters$ret_clu_15 + encounters$ret_clu_30
  encounters$obs_wt_ret_ran <- encounters$total_ret_ran/mean(encounters$total_ret_ran)
  encounters$obs_wt_ret_clu <- encounters$total_ret_clu/mean(encounters$total_ret_clu)

  encounters$cpue_ran_15 <- encounters$ret_ran_15/1400
  encounters$cpue_ran_30 <- encounters$ret_ran_30/2800
  encounters$cpue_clu_15 <- encounters$ret_clu_15/1400
  encounters$cpue_clu_30 <- encounters$ret_clu_30/2800
  
  
  # Compare haul level selectivity ratio method to count models GLM -----
  
  # Poisson models ----
  pois_ran_dat <- 
    encounters |>
    dplyr::select(draw, ret_ran_30, ret_ran_15) |>
    tidyr::pivot_longer(cols = c("ret_ran_30", "ret_ran_15"), values_to = "count") |>
    dplyr::mutate(duration_fac = stringr::str_extract(string = name, pattern = "(\\d)+"),
                  duration_num = as.numeric(duration_fac)*2800/30)
  
  mod_pois <- glm(
    count ~ duration_fac + offset(log(duration_num)), 
    family = poisson(link = "log"), 
    data = pois_ran_dat
  )
  
  pois_clu_dat <- 
    encounters |> 
    dplyr::select(draw, ret_clu_30, ret_clu_15) |>
    tidyr::pivot_longer(cols = c("ret_clu_30", "ret_clu_15"), values_to = "count") |>
    dplyr::mutate(duration_fac = stringr::str_extract(string = name, pattern = "(\\d)+"),
                  duration_num = as.numeric(duration_fac)*2800/30)
  
  mod_pois_clu <- glm(
    count ~ duration_fac + offset(log(duration_num)), 
    family = poisson(link = "log"), 
    data = pois_clu_dat
  )
  
  pois_pred <-
    predict(mod_pois,
            newdata =
              data.frame(duration_fac = c("15", "30"),
                         duration_num = c(1400, 2800)),
            type = "response")
  
  pois_pred_clust <-
    predict(mod_pois_clu,
            newdata =
              data.frame(duration_fac = c("15", "30"),
                         duration_num = c(1400, 2800)),
            type = "response"
    )
  
  results$poisson_ran[jj] <- pois_pred[1]/pois_pred[2]
  
  results$poisson_clu[jj] <- pois_pred_clust[1]/pois_pred_clust[2]
  
  # Binomial models (selectivity ratio-type)
  
  encounters$r_15_30_ran <- encounters$ret_ran_15/1400 / (encounters$ret_ran_30/2800 + encounters$ret_ran_15/1400)
  
  encounters$r_15_30_clu <- encounters$ret_clu_15/1400 / (encounters$ret_clu_30/2800 + encounters$ret_clu_15/1400)
  
  mod_binom_ccr_ran <- glm(
    r_15_30_ran ~ 1, family = binomial(link = "logit"), 
    data = encounters
    )
  
  mod_binom_ccr_clu <- glm(
    r_15_30_clu ~ 1, 
    family = binomial(link = "logit"), 
    data = encounters)
  
  results$binom_ccr_ran[jj] <- predict(mod_binom_ccr_ran, newdata = data.frame(x = 1), type = "response")
  
  results$binom_ccr_clu[jj] <- predict(mod_binom_ccr_clu, newdata = data.frame(x = 1), type = "response")
  
  # Binomial models with observation weighting (selectivity ratio-type)
  
  mod_binom_ccr_wt_ran <- glm(
    r_15_30_ran ~ 1, family = binomial(link = "logit"), 
    weight = obs_wt_ret_ran, 
    data = encounters
  )
  
  mod_binom_ccr_wt_clu <- glm(
    r_15_30_clu ~ 1, 
    family = binomial(link = "logit"), 
    weight = obs_wt_ret_clu, 
    data = encounters
  )
  
  results$binom_ccr_wt_ran[jj] <- predict(mod_binom_ccr_wt_ran, newdata = data.frame(x = 1), type = "response")
  
  results$binom_ccr_wt_clu[jj] <- predict(mod_binom_ccr_wt_clu, newdata = data.frame(x = 1), type = "response")
  
  # Aggregate selectivity ratios ----
  
  results$agg_ccr_ran[jj] <- sum(encounters$ret_ran_15)/1400 / (sum(encounters$ret_ran_30)/2800 + sum(encounters$ret_ran_15)/1400)
  
  results$agg_ccr_clu[jj] <- sum(encounters$ret_clu_15)/1400 / (sum(encounters$ret_clu_30)/2800 + sum(encounters$ret_clu_15)/1400)
  
}

end_time <- Sys.time()

difftime(end_time, start_time)


ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_freqpoly(data = results |>
                  dplyr::filter(cluster_density_n_km2 == 1, cluster_radius_m == 400) |>
                   dplyr::select(fish_density_n_km2, n_hauls, poisson_ran, binom_ccr_ran) |>
                   tidyr::pivot_longer(cols = c("poisson_ran", "binom_ccr_ran")),
                 mapping = aes(x = value, color = name), alpha = 0.5, binwidth = 0.03, position = "dodge") +
  scale_y_continuous(name = "Frequency") +
  scale_x_continuous(name = "Catch ratio", breaks = seq(0, 1, 0.25)) +
  ggtitle(label = "Randomly distributed fish ") +
  facet_grid(n_hauls~fish_density_n_km2) +
  theme_bw()


ggplot() +
  geom_density_ridges(
    data = results |>
      dplyr::filter(cluster_density_n_km2 == 0.5, cluster_radius_m == 1000) |>
      dplyr::select(fish_density_n_km2, n_hauls, poisson_ran, binom_ccr_ran, binom_ccr_wt_ran, poisson_clu, binom_ccr_clu, binom_ccr_wt_clu, agg_ccr_ran, agg_ccr_clu) |>
      tidyr::pivot_longer(cols = c("poisson_ran", "binom_ccr_ran",  "binom_ccr_wt_ran", "poisson_clu", "binom_ccr_clu", "binom_ccr_wt_clu", "agg_ccr_ran", "agg_ccr_clu")),
    mapping = 
      aes(x = value, 
          y = name,
          fill = 
            name), 
    alpha = 0.4, 
    color = "grey40",
    calc_ecdf = TRUE,
    quantiles = 0.5,
    quantile_lines = TRUE
  ) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  scale_y_discrete(name = "Density", limits = rev) +
  scale_x_continuous(name = "Catch comparison rate", breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_fill_viridis_d(option = "H") +
  # ggtitle(label = "Randomly distributed fish ") +
  facet_grid(n_hauls~fish_density_n_km2) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))

ggplot() +
  geom_density_ridges(
    data = results |>
      dplyr::filter(n_hauls == 40) |>
      dplyr::select(fish_density_n_km2, cluster_density_n_km2, cluster_radius_m, n_hauls, binom_ccr_ran, poisson_clu, binom_ccr_clu, binom_ccr_wt_clu, agg_ccr_clu) |>
      tidyr::pivot_longer(cols = c( "poisson_clu", "binom_ccr_clu", "binom_ccr_wt_clu", "agg_ccr_clu")),
    mapping = 
      aes(x = value, 
          y = paste0(name, "-", format(cluster_density_n_km2, nsmall = 1)),
          fill = 
            name), 
    alpha = 0.4, 
    color = "grey40",
    calc_ecdf = TRUE,
    quantiles = 0.5,
    quantile_lines = TRUE
  ) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  scale_y_discrete(name = "Density", limits = rev) +
  scale_x_continuous(name = "Catch comparison rate", breaks = seq(0, 1, 0.25), limits = c(0,1)) +
  scale_fill_viridis_d(option = "H") +
  facet_grid(cluster_radius_m~fish_density_n_km2) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))

ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_density_ridges(
    data = results |>
      dplyr::filter(cluster_density_n_km2 == 2, cluster_radius_m == 400) |>
      dplyr::select(fish_density_n_km2, n_hauls, poisson_ran, binom_ccr_ran, poisson_clu, binom_ccr_clu, agg_ccr_ran, agg_ccr_clu) |>
      tidyr::pivot_longer(cols = c("poisson_ran", "binom_ccr_ran", "poisson_clu", "binom_ccr_clu", "agg_ccr_ran", "agg_ccr_clu")),
    mapping = 
      aes(x = value, 
          y = name,
          fill = 
            name), 
    alpha = 0.6, 
    calc_ecdf = TRUE,
    quantiles = 0.5,
    quantile_lines = TRUE
  ) +
  scale_y_discrete(name = "Sample/Model") +
  scale_x_continuous(name = "Catch comparison rate", breaks = seq(0, 1, 0.25)) +
  scale_fill_viridis_d(option = "H") +
  # ggtitle(label = "Randomly distributed fish ") +
  facet_grid(n_hauls~fish_density_n_km2) + 
  theme_bw()

ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_density_ridges(
    data = results |>
      dplyr::filter(cluster_density_n_km2 == 2, n_hauls  == 80) |>
      dplyr::select(fish_density_n_km2, cluster_radius_m, n_hauls, poisson_ran, binom_ccr_ran, poisson_clu, binom_ccr_clu, agg_ccr_ran, agg_ccr_clu) |>
      tidyr::pivot_longer(cols = c("poisson_ran", "binom_ccr_ran", "poisson_clu", "binom_ccr_clu", "agg_ccr_ran", "agg_ccr_clu")),
    mapping = 
      aes(x = value, 
          y = name,
          fill = 
            name), 
    alpha = 0.6, 
    calc_ecdf = TRUE,
    quantiles = 0.5,
    quantile_lines = TRUE
  ) +
  scale_y_discrete(name = "Sample/Model") +
  scale_x_continuous(name = "Catch comparison rate", breaks = seq(0, 1, 0.25)) +
  scale_fill_viridis_d(option = "H") +
  # ggtitle(label = "Randomly distributed fish ") +
  facet_grid(cluster_radius_m~fish_density_n_km2) + 
  theme_bw()

ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_density_ridges(
    data = results |>
      dplyr::filter(cluster_density_n_km2 == 2, n_hauls == 80, fish_density_n_km2) |>
      dplyr::select(fish_density_n_km2, cluster_radius_m, n_hauls, poisson_ran, binom_ccr_ran, poisson_clu, binom_ccr_clu, agg_ccr_ran, agg_ccr_clu) |>
      tidyr::pivot_longer(cols = c("poisson_ran", "binom_ccr_ran", "poisson_clu", "binom_ccr_clu", "agg_ccr_ran", "agg_ccr_clu")),
    mapping = 
      aes(x = value, 
          y = name,
          fill = 
            name), 
    alpha = 0.6, 
    calc_ecdf = TRUE,
    quantiles = 0.5,
    quantile_lines = TRUE
  ) +
  scale_y_discrete(name = "Sample/Model") +
  scale_x_continuous(name = "Catch comparison rate", breaks = seq(0, 1, 0.25)) +
  scale_fill_viridis_d(option = "H") +
  # ggtitle(label = "Randomly distributed fish ") +
  facet_grid(cluster_radius_m~fish_density_n_km2) + 
  theme_bw()


ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = poisson_clu)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)

ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = binom_ccr_clu)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)

ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = poisson_ran)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)

ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = binom_ccr_ran)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)


  print(ggplot() +
          geom_point(data = poisson_sample$cluster_points[[80]],
                     mapping = aes(x = x_m, y = y_m))) +
    scale_x_continuous(limits = c(3000, 1500))
  
  print(ggplot() +
          geom_point(data = poisson_sample$random_points[[80]],
                     mapping = aes(x = x_m, y = y_m))) +
    scale_x_continuous(limits = c(3000, 1500))

