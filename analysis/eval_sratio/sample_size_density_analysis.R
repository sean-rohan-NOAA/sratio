# Plot effort reduction scenarios
library(sratio)
library(akgfmaps)
library(ggridges)
library(cowplot)


seed <- NULL
ref_distance_fished_m <- 2800
distance_between_tows_m <- 600
# fish_density_n_km2 = c(30, 50, 100, 150, 250, 500, 1000, 2000)
grid_dim_m = c(3000, 1500)
cluster_density_n_km2 = 2
cluster_radius_m = 400
open_boundary = TRUE
# draws = 100

results <- 
  expand.grid(
    sample = 1:1000,
    fish_density_n_km2 = c(25, 50, 100, 200, 500, 1000, 2000),
    n_hauls = c(40, 50, 80), 
    retention = 0.9,
    poisson_ran = NA,
    poisson_clu = NA,
    sratio_ran = NA,
    sratio_clu = NA
  )


start_time <- Sys.time()

for(jj in 1:nrow(results)) {
  
  print(paste0(jj, " - ", Sys.time()))
  
  poisson_sample <- 
    try(sim_matern_clusters(
      fish_density_n_km2 = results$fish_density_n_km2[jj],
      grid_dim_m = grid_dim_m,
      cluster_density_n_km2 = cluster_density_n_km2,
      cluster_radius_m = cluster_radius_m,
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
  
  # Simulate retention
  set.seed(seed)
  
  for(ii in 2:ncol(encounters)) {
    
    encounters <- cbind(
      encounters, 
      mapply(function(n) rbinom(1, n, results$retention[jj]), encounters[, ii])
    )
    
    names(encounters)[ncol(encounters)] <- gsub(x = names(encounters)[ii], pattern = "enc", replacement = "ret")
    
  }
  
  
  encounters$cpue_ran_15 <- encounters$ret_ran_15/1400
  encounters$cpue_ran_30 <- encounters$ret_ran_30/2800
  encounters$cpue_clu_15 <- encounters$ret_clu_15/1400
  encounters$cpue_clu_30 <- encounters$ret_clu_30/2800
  
  test_dat <- 
    encounters |>
    dplyr::mutate(
      sample = results$sample[jj],
      n_hauls = results$n_hauls[jj],
      fish_density_n_km2 = results$fish_density_n_km2[jj],
      retention = results$retention[jj]
    ) |>
    dplyr::group_by(sample, n_hauls, fish_density_n_km2, retention) |>
    dplyr::summarise(
      total_enc_ran_30 = sum(enc_ran_30, na.rm = TRUE),
      total_enc_ran_15 = sum(enc_ran_15, na.rm = TRUE), 
      total_enc_clu_30 = sum(enc_clu_30, na.rm = TRUE), 
      total_enc_clu_15 = sum(enc_clu_15, na.rm = TRUE), 
      total_ret_ran_30 = sum(ret_ran_30, na.rm = TRUE), 
      total_ret_ran_15 = sum(ret_ran_15, na.rm = TRUE), 
      total_ret_clu_30 = sum(ret_clu_30, na.rm = TRUE), 
      total_ret_clu_15 = sum(ret_clu_15, na.rm = TRUE),
      .groups = "keep"
    ) |>
    dplyr::mutate(effort_30 = 2800 * results$n_hauls[jj],
                  effort_15 = 1400 * results$n_hauls[jj])
  
  
  
  # Compare haul level selectivity ratio method to count models GLM
  
  # Poisson models ----
  mod_pois <- glm(value ~ duration_fac + offset(log(duration_num)), 
                  family = poisson(link = "log"), 
                  data = 
                    dplyr::select(encounters, draw, ret_ran_30, ret_ran_15) |>
                    tidyr::pivot_longer(cols = c("ret_ran_30", "ret_ran_15")) |>
                    dplyr::mutate(duration_fac = stringr::str_extract(string = name, pattern = "(\\d)+"),
                                  duration_num = as.numeric(duration_fac)*2800/30)
  )
  
  mod_pois_clu <- glm(value ~ duration_fac + offset(log(duration_num)), 
                        family = poisson(link = "log"), 
                        data =  
                          dplyr::select(encounters, draw, ret_clu_30, ret_clu_15) |>
                          tidyr::pivot_longer(cols = c("ret_clu_30", "ret_clu_15")) |>
                          dplyr::mutate(duration_fac = stringr::str_extract(string = name, pattern = "(\\d)+"),
                                        duration_num = as.numeric(duration_fac)*2800/30))
  
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
            type = "response")
  
  results$poisson_ran[jj] <- pois_pred[1]/pois_pred[2]
  
  results$poisson_clu[jj] <- pois_pred_clust[1]/pois_pred_clust[2]
  
  # Binomial sratio models
  encounters$r_15_30 <- encounters$enc_ran_15/1400 / (encounters$enc_ran_30/2800 + encounters$enc_ran_15/1400)
  
  encounters$r_15_30_clu <- encounters$enc_clu_15/1400 / (encounters$enc_clu_30/2800 + encounters$enc_clu_15/1400)
  
  mod_sratio <- glm(r_15_30 ~ 1, family = binomial(link = "logit"), data = encounters)
  
  mod_sratio_clu <- glm(r_15_30_clu ~ 1, family = binomial(link = "logit"), data = encounters)northeno
  
  results$sratio_ran[jj] <- predict(mod_sratio, newdata = data.frame(x = 1), type = "response")
  
  results$sratio_clu[jj] <- predict(mod_sratio_clu, newdata = data.frame(x = 1), type = "response")
  
}

end_time <- Sys.time()

difftime(end_time, start_time)


ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_freqpoly(data = results |>
                   dplyr::select(fish_density_n_km2, n_hauls, poisson_ran, sratio_ran) |>
                   tidyr::pivot_longer(cols = c("poisson_ran", "sratio_ran")),
                 mapping = aes(x = value, color = name), alpha = 0.5, binwidth = 0.03, position = "dodge") +
  scale_y_continuous(name = "Frequency") +
  scale_x_continuous(name = "Catch ratio", breaks = seq(0, 1, 0.25)) +
  ggtitle(label = "Randomly distributed fish ") +
  facet_grid(n_hauls~fish_density_n_km2) + 
  theme_bw()


ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_density_ridges(
    data = results |>
      dplyr::select(fish_density_n_km2, n_hauls, poisson_ran, sratio_ran, poisson_clu, sratio_clu) |>
      tidyr::pivot_longer(cols = c("poisson_ran", "sratio_ran", "poisson_clu", "sratio_clu")),
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
  scale_x_continuous(name = "Catch ratio", breaks = seq(0, 1, 0.25)) +
  scale_fill_viridis_d(option = "H") +
  # ggtitle(label = "Randomly distributed fish ") +
  facet_grid(n_hauls~fish_density_n_km2) + 
  theme_bw()


ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = poisson_clu)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)

ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = sratio_clu)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)

ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = poisson_ran)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)

ggplot() +
  geom_histogram(data = results,
                 mapping = aes(x = sratio_ran)) +
  geom_vline(xintercept = 0.5, linetype = 2) +
  facet_grid(n_hauls~fish_density_n_km2)
