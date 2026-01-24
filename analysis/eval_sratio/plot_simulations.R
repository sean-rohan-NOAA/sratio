library(cowplot)

ex_analysis <- 
  expand.grid(
    fish_density_n_km2 = c(20, 50, 100, 500, 1000),
    cluster_density = c(0.5, 1, 2),
    cluster_radius_m = c(400, 1000, 2000)
  )

random_points <- data.frame()
cluster_points <- data.frame()

for(ii in 1:nrow(ex_analysis)) {
  
  sim_dist <- sim_matern_clusters(
    fish_density_n_km2 = ex_analysis$fish_density_n_km2[ii],
    grid_dim_m = c(3000, 3000),
    cluster_density_n_km2 = ex_analysis$cluster_density[ii],
    cluster_radius_m = ex_analysis$cluster_radius_m[ii],
    open_boundary = TRUE,
    draws = 1
  )
  
  ran <- as.data.frame(sim_dist$random_points[[1]])
  clu <- as.data.frame(sim_dist$cluster_points[[1]])
  
  ran$type <- "Random"
  ran$fish_density_n_km2 <- ex_analysis$fish_density_n_km2[ii]

  
  clu$type <- "Cluster"
  clu$fish_density_n_km2 <- ex_analysis$fish_density_n_km2[ii]
  clu$cluster_density_n_km2 <- ex_analysis$cluster_density[ii]
  clu$cluster_radius_m <- ex_analysis$cluster_radius_m[ii]

  if(ex_analysis$cluster_density[ii] == 0.5 & ex_analysis$cluster_radius_m[ii] == 400) {
    random_points <- rbind(random_points, ran)
  }

  cluster_points <- rbind(cluster_points, clu)
  
}

random_points <- unique(random_points)

ex_sim_points <- dplyr::bind_rows(random_points, cluster_points)


dplyr::select(random_points, type, fish_density_n_km2, draw) |>
  unique()

ggplot() +
  geom_point(data = 
               random_points,
             mapping = aes(x = x_m/1000, y = y_m/1000),
             size = rel(0.2),
             color = "grey70") +
  geom_polygon(data = data.frame(x = c(0, 2.8, 2.8, 0, 0), y = c(0.1, 0.1, .117, .117, .1)),
               mapping = aes(x = x, y = y), fill = "red", alpha = 0.5) +
  geom_polygon(data = data.frame(x = c(0, 1.4, 1.4, 0, 0), y = c(0.7, 0.7, 0.717, 0.717, 0.7)),
               mapping = aes(x = x, y = y), fill = "blue", alpha = 0.5) +
  facet_wrap(~fish_density_n_km2) +
  scale_x_continuous(name = "X (km)", limits = c(0, 3), expand = c(0, 0)) +
  scale_y_continuous(name = "Y (km)", limits = c(0, 3), expand = c(0, 0)) +
  theme_few() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggplot() +
  geom_point(data = 
               dplyr::filter(cluster_points, cluster_density_n_km2 == 0.5),
             mapping = aes(x = x_m/1000, y = y_m/1000),
             size = rel(0.2)) +
  geom_polygon(data = data.frame(x = c(0, 2.8, 2.8, 0, 0), y = c(0.1, 0.1, 0.117, 0.117, 0.1)),
               mapping = aes(x = x, y = y), fill = "red", alpha = 0.5) +
  geom_polygon(data = data.frame(x = c(0, 1.4, 1.4, 0, 0), y = c(0.7, 0.7, 0.717, 0.717, 0.7)),
               mapping = aes(x = x, y = y), fill = "blue", alpha = 0.5) +
  facet_grid(cluster_radius_m~fish_density_n_km2) +
  scale_x_continuous(name = "X (km)", limits = c(0, 3), expand = c(0, 0)) +
  scale_y_continuous(name = "Y (km)", limits = c(0, 3), expand = c(0, 0)) +
  theme_few() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggplot() +
  geom_point(data = 
               dplyr::filter(cluster_points, fish_density_n_km2 == 500),
             mapping = aes(x = x_m/1000, y = y_m/1000),
             size = rel(0.2)) +
  geom_polygon(data = data.frame(x = c(0, 2.8, 2.8, 0, 0), y = c(0.1, 0.1, 0.117, 0.117, 0.1)),
               mapping = aes(x = x, y = y), fill = "red", alpha = 0.5) +
  geom_polygon(data = data.frame(x = c(0, 1.4, 1.4, 0, 0), y = c(0.7, 0.7, 0.717, 0.717, 0.7)),
               mapping = aes(x = x, y = y), fill = "blue", alpha = 0.5) +
  facet_grid(cluster_radius_m~cluster_density_n_km2) +
  scale_x_continuous(name = "X (km)", limits = c(0, 3), expand = c(0, 0)) +
  scale_y_continuous(name = "Y (km)", limits = c(0, 3), expand = c(0, 0)) +
  theme_few() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))






