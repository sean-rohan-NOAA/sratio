
seed <- NULL
ref_distance_fished_m <- 2800
distance_between_tows_m <- 1000
fish_density_n_km2 = 150
grid_dim_m = c(3000, 3000)
cluster_density_n_km2 = 2
cluster_radius_m = 400
open_boundary = TRUE
draws = 100
n_iter = 1000

# Function to generate distributions
sim_matern_clusters <- 
  function(
    fish_density_n_km2,
    grid_dim_m,
    cluster_density_n_km2,
    cluster_radius_m,
    open_boundary = TRUE,
    seed, 
    draws
  ) {
    
    cluster_points <- vector(mode = "list", length = draws)
    random_points <- vector(mode = "list", length = draws)
    
    set.seed(seed)
    
    for(ii in 1:draws) {
      
      # Simulate Poisson point process ----
      
      # Poisson fish
      n_points <- rpois(1, fish_density_n_km2 * prod(grid_dim_m) / 1e6)
      
      # Assign coordinates to points
      x_m <- grid_dim_m[1] * runif(n_points)
      
      y_m <- grid_dim_m[2] * runif(n_points)
      
      draw <- rep(ii, n_points)
      
      random_points[[ii]] <-
        cbind(x_m,
              y_m,
              draw)
      
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
      x_m <- rep(x_cluster, cluster_sizes) + xx0
      
      y_m <- rep(y_cluster, cluster_sizes) + yy0
      
      cluster_points[[ii]] <- 
        cbind(x_m,
              y_m,
              draw)
      
    }
    
    output <- 
      list(random_points = random_points,
           cluster_points = cluster_points)
    
    return(output)
    
  }


# Function to count count encounters
sim_encounters <- 
  function(x,
           origin_m,
           effort_m,
           make_plot = FALSE) {
    
    # Points that overlap with trawl path
    encounters <- sum(
      (x[,1] >= origin_m[1] & x[,1] <= origin_m[1] + effort_m[1]) & 
        (x[,2] >= origin_m[2] & x[,2] <= origin_m[2] + effort_m[2])
    )
    
    return(encounters)
  
  }


for(jj in 1:length(n_iter)) {
  
  test <- 
    sim_matern_clusters(
      fish_density_n_km2 = fish_density_n_km2,
      grid_dim_m = grid_dim_m,
      cluster_density_n_km2 = cluster_density_n_km2,
      cluster_radius_m = cluster_radius_m,
      open_boundary = open_boundary,
      seed = seed, 
      draws = draws
    )
  
  
  encounters <-
    data.frame(
      draw = 1:draws,
      enc_ran_30 = 
        unlist(
          lapply(
            X = test[['random_points']], 
            FUN = sim_encounters, 
            origin_m = c(0, grid_dim_m[2]/2), 
            effort_m = c(ref_distance_fished_m, 17)
          )
        ),
      enc_ran_15 = 
        unlist(
          lapply(
            X = test[['random_points']], 
            FUN = sim_encounters, 
            origin_m = c(0, grid_dim_m[2]/2 + distance_between_tows_m), 
            effort_m = c(ref_distance_fished_m/2, 17)
          )
        ),
      enc_ran_5 = 
        unlist(
          lapply(
            X = test[['random_points']], 
            FUN = sim_encounters, 
            origin_m = c(0, grid_dim_m[2]/2 + distance_between_tows_m), 
            effort_m = c(ref_distance_fished_m/6, 17)
          )
        ),
      enc_clu_30 = 
        unlist(
          lapply(
            X = test[['cluster_points']], 
            FUN = sim_encounters, 
            origin_m = c(0, grid_dim_m[2]/2), 
            effort_m = c(ref_distance_fished_m, 17)
          )
        ),
      enc_clu_15 = 
        unlist(
          lapply(
            X = test[['cluster_points']], 
            FUN = sim_encounters, 
            origin_m = c(0, grid_dim_m[2]/2 + distance_between_tows_m), 
            effort_m = c(ref_distance_fished_m/2, 17)
          )
        ),
      enc_clu_5 = 
        unlist(
          lapply(
            X = test[['cluster_points']], 
            FUN = sim_encounters, 
            origin_m = c(0, grid_dim_m[2]/2 + distance_between_tows_m), 
            effort_m = c(ref_distance_fished_m/6, 17)
          )
        )
    )
  
  
  # Simulate retention
  set.seed(seed)
  
  for(ii in 2:ncol(encounters)) {
    
    encounters <- cbind(
      encounters, 
      mapply(function(n) rbinom(1, n, 0.5), encounters[, ii])
    )
    
    names(encounters)[ncol(encounters)] <- gsub(x = names(encounters)[ii], pattern = "enc", replacement = "ret")
    
  }
  
  
  # Compare selectivity ratio approach to Poisson GLM
  ex_dat <- dplyr::select(encounters, draw, ret_ran_30, ret_ran_15) |>
    tidyr::pivot_longer(cols = c("ret_ran_30", "ret_ran_15")) |>
    dplyr::mutate(duration_fac = stringr::str_extract(string = name, pattern = "(\\d)+"),
                  duration_num = as.numeric(duration_fac)*2800/30)
  
  
  # Pred Poisson
  mod_pois <- glm(value ~ duration_fac + offset(log(duration_num)), family = poisson(link = "log"), data = ex_dat)
  
  pois_pred <-
    predict(mod_pois,
            newdata =
              data.frame(duration_fac = c("15", "30"),
                         duration_num = c(1400, 2800)),
            type = "response")
  
  # Pred sratio binomial
  encounters$r_15_30 <- encounters$enc_ran_15/1400 / (encounters$enc_ran_30/2800 + encounters$enc_ran_15/1400)
  
  mod_sratio <- glm(r_15_30 ~ 1, family = binomial(link = "logit"), data = encounters)
  
  
  
  
  poisson_results[jj] <- pois_pred[1]/pois_pred[2]
  sratio_results[jj] <- predict(mod_sratio, newdata = data.frame(x = 1), type = "response")
  
  
  
}

par(mfrow = c(1, 2))
hist(poisson_results)
abline(v = 0.5, col = "red")

hist(sratio_results)
abline(v = 0.5, col = "red")
