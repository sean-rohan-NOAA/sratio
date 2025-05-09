# Sean Rohan
library(ggplot2)
library(tidyr)
library(ggthemes)

# Function to generate distributions
sim_matern_clusters <- 
  function(
    fish_density_n_km2,
    grid_dim_m,
    cluster_density_n_km2,
    cluster_radius_m,
    open_boundary = TRUE,
    # origin_m_1,
    # origin_m_2,
    # effort_m_1,
    # effort_m_2,
    seed, 
    draws
    # make_plot
  ) {
    
    cluster_points <- vector(mode = "list", length = draws)
    random_points <- vector(mode = "list", length = draws)
    
    set.seed(seed)
    
    for(ii in 1:draws) {
      
      # Simulate Poisson point process ----
      
      # Poisson fish
      n_points <- rpois(1, fish_density_n_km2 * prod(grid_dim_m) / 1e6)
      
      # Assign coordinates to points
      xx_ran <- grid_dim_m[1] * runif(n_points)
      
      yy_ran <- grid_dim_m[2] * runif(n_points)
      
      # Simulate Matern cluster Poisson point process ----
      
      # Poisson clusters
      n_clusters <- rpois(1, cluster_density_n_km2 * prod(grid_dim_m) / 1e6)
      
      # Coordinates for cluster centers
      if(open_boundary) {
        # Case where points can 'spill' over the boundary
        x_cluster <- grid_dim_m[1] * runif(n_clusters)
        y_cluster <- grid_dim_m[2] * runif(n_clusters)
      } else {
        # Case where clusters are seeded at least the distance of the radius away from the boundary
        x_cluster <- cluster_radius_m + (grid_dim_m[1] - 2 * cluster_radius_m) * runif(n_clusters)
        y_cluster <- cluster_radius_m + (grid_dim_m[2] - 2 * cluster_radius_m) * runif(n_clusters)
      }

      # Randomly allocate points to clusters
      cluster_sizes <- tabulate(sample(x = 1:n_clusters, size = n_points, replace = TRUE))
      
      # Generate relative polar coordinates (distance and angle) for each point
      theta <- 2 * pi * runif(n_points)
      
      rho <- cluster_radius_m * sqrt(runif(n_points))
      
      # Convert polar coordinates to Cartesian coordinates
      xx0 <- rho*cos(theta)
      
      yy0 <- rho*sin(theta)
      
      # Translate relative point coordinates to clusters centers
      xx_clus <- rep(x_cluster, cluster_sizes) + xx0
      
      yy_clus <- rep(y_cluster, cluster_sizes) + yy0
      
      
      
      data.frame(xx_ran,
                 yy_ran)
      
      data.frame(xx_clus,
                 yy_clus)
      
      # encounter_random_1  <- 
      #   (xx_ran >= origin_m_1[1] & xx_ran <= origin_m_1[1] + effort_m_1[1]) & 
      #   (yy_ran >= origin_m_1[2] & yy_ran <= origin_m_1[2] + effort_m_1[2])
      # 
      # encounter_random_2  <- 
      #   (xx_ran >= origin_m_2[1] & xx_ran <= origin_m_2[1] + effort_m_2[1]) & 
      #   (yy_ran >= origin_m_2[2] & yy_ran <= origin_m_2[2] + effort_m_2[2])
      # 
      # encounter_clustered_1  <- 
      #   (xx_clus >= origin_m_1[1] & xx_clus <= origin_m_1[1] + effort_m_1[1]) & 
      #   (yy_clus >= origin_m_1[2] & yy_clus <= origin_m_1[2] + effort_m_1[2])
      # 
      # encounter_clustered_2  <- 
      #   (xx_clus >= origin_m_2[1] & xx_clus <= origin_m_2[1] + effort_m_2[1]) & 
      #   (yy_clus >= origin_m_2[2] & yy_clus <= origin_m_2[2] + effort_m_2[2])
      # 
      # encounters[[ii]] <-
      #   data.frame(
      #     draw = ii,
      #     distribution_type = c(rep("random", 2), rep("clustered", 2)),
      #     set_grid_dim_m_x = grid_dim_m[1],
      #     set_grid_dim_m_y = grid_dim_m[2],
      #     set_fish_density_n_km2 = fish_density_n_km2,
      #     set_cluster_density_n_km2 = c(rep(NA, 2), rep(cluster_density_n_km2, 2)),
      #     set_cluster_radius_m = c(rep(NA, 2), rep(cluster_radius_m, 2)),
      #     set_effort_x_m = rep(c(effort_m_1[1], effort_m_2[1]), 2),
      #     set_effort_y_m = rep(c(effort_m_1[2], effort_m_2[2]), 2), 
      #     set_effort_m2 = rep(c(prod(effort_m_1), prod(effort_m_2)), 2),
      #     set_open_boundary = open_boundary,
      #     n_fish = n_points,
      #     n_clusters = c(rep(NA, 2), rep(n_clusters, 2)),
      #     fish_density_n_km2 = n_points / rep(c(prod(effort_m_1), prod(effort_m_2)), 2),
      #     n_encounter = 
      #       c(
      #         sum(encounter_random_1), 
      #         sum(encounter_random_2), 
      #         sum(encounter_clustered_1), 
      #         sum(encounter_clustered_2)
      #       )
      #   )
      
    }
    
    # if(make_plot) {
    #   
    #   # Plot random and clustered results with effort overlaid
    #   effort_m_1_bbox <- 
    #     cbind(
    #       rep(origin_m_1[1], 5) + c(0, effort_m_1[1], effort_m_1[1], 0, 0),
    #       rep(origin_m_1[2], 5) + c(0, 0, effort_m_1[2], effort_m_1[2], 0)
    #     )
    #   
    #   effort_m_2_bbox <-
    #     cbind(
    #       rep(origin_m_2[1], 5) + c(0, effort_m_2[1], effort_m_2[1], 0, 0),
    #       rep(origin_m_2[2], 5) + c(0, 0, effort_m_2[2], effort_m_2[2], 0)
    #     )
    #   
    #   par(mfrow = c(1,2))
    #   plot(
    #     xx_ran, 
    #     yy_ran, 
    #     xlim = c(0, grid_dim_m[1]), 
    #     ylim = c(0, grid_dim_m[2]),
    #     cex = 0.2,
    #     xlab = NA,
    #     ylab = NA,
    #     main = "Random"
    #   )
    #   polygon(x = effort_m_1_bbox[,1], effort_m_1_bbox[,2], col = "red", border = NA)
    #   polygon(x = effort_m_2_bbox[,1], effort_m_2_bbox[,2], col = "blue", border = NA)
    #   
    #   
    #   plot(
    #     xx_clus, 
    #     yy_clus, 
    #     xlim = c(0, grid_dim_m[1]), 
    #     ylim = c(0, grid_dim_m[2]),
    #     cex = 0.2,
    #     xlab = NA,
    #     ylab = NA,
    #     main = "Clustered"
    #   )
    #   polygon(x = effort_m_1_bbox[,1], effort_m_1_bbox[,2], col = "red", border = NA)
    #   polygon(x = effort_m_2_bbox[,1], effort_m_2_bbox[,2], col = "blue", border = NA)
    #   
    # }
    
    output <- do.call(rbind, encounters)
    
    return(output)
    
  }

sim_encounters <-
  function() {
    # origin_m_1,
    # origin_m_2,
    # effort_m_1,
    # effort_m_2,
    # make_plot
  }



encounter_random_1  <- 
  (xx_ran >= origin_m_1[1] & xx_ran <= origin_m_1[1] + effort_m_1[1]) & 
  (yy_ran >= origin_m_1[2] & yy_ran <= origin_m_1[2] + effort_m_1[2])

encounter_random_2  <- 
  (xx_ran >= origin_m_2[1] & xx_ran <= origin_m_2[1] + effort_m_2[1]) & 
  (yy_ran >= origin_m_2[2] & yy_ran <= origin_m_2[2] + effort_m_2[2])

encounter_clustered_1  <- 
  (xx_clus >= origin_m_1[1] & xx_clus <= origin_m_1[1] + effort_m_1[1]) & 
  (yy_clus >= origin_m_1[2] & yy_clus <= origin_m_1[2] + effort_m_1[2])

encounter_clustered_2  <- 
  (xx_clus >= origin_m_2[1] & xx_clus <= origin_m_2[1] + effort_m_2[1]) & 
  (yy_clus >= origin_m_2[2] & yy_clus <= origin_m_2[2] + effort_m_2[2])

encounters[[ii]] <-
  data.frame(
    draw = ii,
    distribution_type = c(rep("random", 2), rep("clustered", 2)),
    set_grid_dim_m_x = grid_dim_m[1],
    set_grid_dim_m_y = grid_dim_m[2],
    set_fish_density_n_km2 = fish_density_n_km2,
    set_cluster_density_n_km2 = c(rep(NA, 2), rep(cluster_density_n_km2, 2)),
    set_cluster_radius_m = c(rep(NA, 2), rep(cluster_radius_m, 2)),
    set_effort_x_m = rep(c(effort_m_1[1], effort_m_2[1]), 2),
    set_effort_y_m = rep(c(effort_m_1[2], effort_m_2[2]), 2), 
    set_effort_m2 = rep(c(prod(effort_m_1), prod(effort_m_2)), 2),
    set_open_boundary = open_boundary,
    n_fish = n_points,
    n_clusters = c(rep(NA, 2), rep(n_clusters, 2)),
    fish_density_n_km2 = n_points / rep(c(prod(effort_m_1), prod(effort_m_2)), 2),
    n_encounter = 
      c(
        sum(encounter_random_1), 
        sum(encounter_random_2), 
        sum(encounter_clustered_1), 
        sum(encounter_clustered_2)
      )
  )


test_densities <- c(50, seq(100, 2000, 100))

eval_density <- vector(mode = "list", length = length(test_densities))

for(jj in 1:length(test_densities)) {
  
  eval_density[[jj]] <-
    sim_encounters_matern_cluster(
      fish_density_n_km2 = test_densities[jj],
      grid_dim_m = c(3000, 3000),
      cluster_density_n_km2 = 2,
      open_boundary = TRUE,
      cluster_radius_m = 400,
      origin_m_1 = c(0, 1000),
      origin_m_2 = c(0, 2000),
      effort_m_1 = c(2800, 17),
      effort_m_2 = c(1400, 17),
      seed = 1337, 
      draws = 200,
      make_plot = TRUE
    )
  
}

results <- 
  do.call(rbind, eval_density) |>
  dplyr::mutate(tow_duration = ifelse(set_effort_x_m == 2800, "long", "short")) |>
  dplyr::select(draw, distribution_type, tow_duration, set_fish_density_n_km2, n_fish, n_encounter) |>
  tidyr::pivot_wider(names_from = tow_duration, values_from = n_encounter) |>
  dplyr::mutate(p = short/long,
                r = short/(long+short))

ggplot() +
  geom_boxplot(data = results, 
               mapping = aes(x = set_fish_density_n_km2,
                             y = p,
                             group = set_fish_density_n_km2)) +
  geom_hline(yintercept = 0.5) +
  facet_wrap(~distribution_type) +
  scale_x_continuous(name = expression('Density (#'%.%'km'^2*')')) +
  scale_y_continuous(name = expression('Encounter ratio, '*n[15]/n[30]), limits = c(0, 5)) +
  # scale_color_colorblind(name = "Distribution") +
  theme_bw()


ggplot() +
  geom_boxplot(data = results, 
               mapping = aes(x = set_fish_density_n_km2,
                             y = r,
                             group = set_fish_density_n_km2)) +
  geom_hline(yintercept = 1/3) +
  facet_wrap(~distribution_type) +
  scale_x_continuous(name = expression('Density (#'%.%'km'^2*')')) +
  scale_y_continuous(name = expression('Encounter comparison rate, '*n[15]/(n[30]+n[15]))) +
  # scale_color_colorblind(name = "Distribution") +
  theme_bw()


