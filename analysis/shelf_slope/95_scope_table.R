library(sratio)

# 2023 Scope Values
scope_table_2023 <- readRDS(here::here("analysis", "shelf_slope", "data", "ss_haul.rds")) |> 
  dplyr::select(dplyr::all_of(c("BOTTOM_DEPTH", "GEAR", "WIRE_LENGTH"))) |>
  dplyr::inner_join(data.frame(GEAR = c(44, 172), 
                               GEAR_NAME = c("83-112", "PNE"), 
                               table = c("EBS shelf", "EBS slope")))

# EBS shelf and slope standard scope tables
scope_table_standard_m <- readxl::read_xlsx(
  here::here("analysis", "shelf_slope", "docs", "shelf_slope_table.xlsx")
  ) |>
  dplyr::select(dplyr::all_of(c("min_depth_m", "max_depth_m", "wire_out_m", "table"))) |>
  tidyr::pivot_longer(cols = "wire_out_m") |>
  dplyr::mutate(mean_depth_m = (min_depth_m+max_depth_m)/2)

scope_table_standard_fm <- readxl::read_xlsx(
  here::here("analysis", "shelf_slope", "docs", "shelf_slope_table.xlsx")
) |>
  dplyr::select(dplyr::all_of(c("min_depth_fm", "max_depth_fm", "wire_out_fm", "table"))) |>
  tidyr::pivot_longer(cols = "wire_out_fm") |>
  dplyr::mutate(mean_depth_fm = (min_depth_fm+max_depth_fm)/2)


# Plot 2023 scope table
ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "scope_tables", "2023_depth_vs_scope.png"),
              width = 6, height = 4, units = "in", res = 180)
print(
  ggplot() +
    geom_point(data = scope_table_2023,
               mapping = aes(x = BOTTOM_DEPTH, 
                             y = WIRE_LENGTH, 
                             color = table, 
                             shape = table)) +
    scale_x_continuous(name = "Bottom Depth (m)") +
    scale_y_continuous(name = "Wire Length (m)") +
    ggtitle(label = "2023 Shelf/Slope Scope") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
)
dev.off()


ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "scope_tables", "2023_depth_vs_scope_w_table.png"),
              width = 6, height = 4, units = "in", res = 180)
print(
  ggplot() +
    geom_point(data = scope_table_2023,
               mapping = aes(x = BOTTOM_DEPTH, 
                             y = WIRE_LENGTH, 
                             color = table, 
                             shape = table),
               size = rel(3)) +
    geom_path(data = scope_table_standard_m,
              mapping = aes(x = mean_depth_m, y = value, color = table), linewidth = 1.2) +
    scale_x_continuous(name = "Bottom Depth (m)") +
    scale_y_continuous(name = "Wire Length (m)") +
    ggtitle(label = "2023 Shelf/Slope Scope") +
    scale_shape(name = "Obs.", solid = FALSE) +
    scale_color_colorblind(name = "Table") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
)
dev.off()


print(
  ggplot() +
    geom_path(data = scope_table_standard_fm,
              mapping = aes(x = mean_depth_fm, y = value, color = table), linewidth = 1.2) +
    scale_x_continuous(name = "Bottom Depth (m)") +
    scale_y_continuous(name = "Wire Length (m)") +
    ggtitle(label = "2023 Shelf/Slope Scope") +
    scale_shape(name = "Obs.", solid = FALSE) +
    scale_color_colorblind(name = "Table") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
)


# Create new scope tables

estimate_new_scope <- function(target_scope, scope, min_z, max_z) {
  
  index_upr <- min(which(target_scope < scope))
  index_lwr <- max(which(target_scope > scope))

  diff_scope <- scope[index_upr] - scope[index_lwr]
  mean_z_upr <- (min_z[index_upr] + max_z[index_upr]-1)/2
  mean_z_lwr <- (min_z[index_lwr] + max_z[index_lwr]+1)/2
  
  diff_z <- mean_z_upr - mean_z_lwr
  
             
  new_mid_z <- mean_z_lwr + diff_z/diff_scope * (target_scope - scope[index_lwr])
  
  new_min_z <- ceiling(new_mid_z - (diff_z/2))
  new_max_z <- floor(new_mid_z + (diff_z/2))
  
  out <- data.frame(scope = target_scope,
                    min_z = new_min_z,
                    max_z = new_max_z,
                    mean_z = (new_min_z + new_max_z)/2)
  
  return(out)
  
}


# Create modified slope scope table
scope_vals_slope_z <- seq(325, 600, 25)

modified_slope_scope <- data.frame(scope = 300,
                                   min_z = 80,
                                   max_z = 93) |>
  dplyr::mutate(mean_z = (min_z+max_z)/2) |>
  dplyr::bind_rows(dplyr::filter(scope_table_standard_fm, 
                                 table == "EBS shelf", 
                                 max_depth_fm < 80) |>
                     dplyr::mutate(
                       max_depth_fm = dplyr::if_else(max_depth_fm == 78, 79, max_depth_fm)
                       ) |>
                     dplyr::rename(min_z = min_depth_fm,
                                   max_z = max_depth_fm,
                                   mean_z = mean_depth_fm,
                                   scope = value)) |>
  dplyr::select(-table, -name) |>
  dplyr::arrange(scope)


for(ii in 1:length(scope_vals_slope_z)) {
  modified_slope_scope <- dplyr::bind_rows(
    modified_slope_scope,
    estimate_new_scope(target_scope = scope_vals_slope_z[ii], 
                       scope = scope_table_standard_fm$value[scope_table_standard_fm$table == "EBS slope"], 
                       min_z = scope_table_standard_fm$min_depth_fm[scope_table_standard_fm$table == "EBS slope"], 
                       max_z = scope_table_standard_fm$max_depth_fm[scope_table_standard_fm$table == "EBS slope"])
  )

}

modified_slope_scope <- modified_slope_scope |>
  dplyr::mutate() |>
  dplyr::rename(min_depth_fm = min_z,
                max_depth_fm = max_z,
                mean_depth_fm = mean_z,
                scope_fm = scope) |>
  dplyr::mutate(table = "EBS slope-modified",
                scope_m = 1.8288 * scope_fm,
                min_depth_m = 1.8288 * min_depth_fm, 
                max_depth_m = 1.8288 * max_depth_fm)


# Create modified shelf scope table
modified_shelf_scope <- 
  dplyr::filter(scope_table_standard_fm, 
                table == "EBS shelf") |>
  dplyr::select(-table, -name, scope_fm = value) |>
  dplyr::bind_rows(
    dplyr::filter(modified_slope_scope, 
                  scope_fm > 350)
  ) |>
  dplyr::mutate(table = "EBS shelf-extended",
                scope_m = 1.8288 * scope_fm,
                min_depth_m = 1.8288 * min_depth_fm,
                max_depth_m = 1.8288 * max_depth_fm
                )

table_order <- c("EBS slope",
                "EBS shelf",
                "GOA/AI",
                "EBS slope-modified", 
                "EBS shelf-extended")


plot_new_table <- ggplot() +
  geom_path(data = scope_table_standard_fm |>
              dplyr::rename(scope_fm = value),
            mapping = aes(x = mean_depth_fm, 
                          y = scope_fm, 
                          color = factor(table,
                                         levels = table_order)),
            linewidth = 1,
            linetype = 3) +
  geom_path(data = dplyr::bind_rows(modified_shelf_scope,
                                    modified_slope_scope) ,
            mapping = aes(x = mean_depth_fm, 
                          y = scope_fm, 
                          color = factor(table,
                                         levels = table_order)),
            linewidth = 1.3,
            linetype = 2) +
  scale_x_continuous(name = "Bottom Depth (fm)") +
  scale_y_continuous(name = "Wire Length (fm)") +
  scale_shape(name = "Obs.", solid = FALSE) +
  scale_color_manual(values = c("EBS slope" = "orange",
                                "EBS shelf" = "#CC79A7",
                                "GOA/AI" = "#009E73",
                                "EBS slope-modified" = "#000000", 
                                "EBS shelf-extended" = "#56B4E9")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.2))

ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "scope_tables", "2024_modified_scope.png"), width = 6, height = 4, units = "in", res = 180)
print(plot_new_table)
dev.off()


dplyr::select(modified_shelf_scope, 
              dplyr::all_of(c("min_depth_fm", "max_depth_fm", "scope_fm", 
                              "min_depth_m", "max_depth_m", "scope_m"))) |>
  write.csv(here::here("analysis", "shelf_slope", "plots", "scope_tables", "2024_modified_scope_shelf_unedited.csv"),
            row.names = FALSE)

dplyr::select(modified_slope_scope, 
              dplyr::all_of(c("min_depth_fm", "max_depth_fm", "scope_fm", 
                              "min_depth_m", "max_depth_m", "scope_m"))) |>
  write.csv(here::here("analysis", "shelf_slope", "plots", "scope_tables", "2024_modified_scope_slope_unedited.csv"),
            row.names = FALSE)

# Check edited scope table

modified_slope_edited <- read.csv(here::here("analysis", "shelf_slope", "plots", "scope_tables", "2024_modified_scope_slope_edited.csv"))

modified_shelf_edited <- read.csv(here::here("analysis", "shelf_slope", "plots", "scope_tables", "2024_modified_scope_shelf_edited.csv"))


plot_slope_table_fm <- ggplot() +
  geom_segment(data = 
                 dplyr::bind_rows(modified_slope_edited |>
                                    dplyr::mutate(table = "EBS slope-modified"),
                                       modified_shelf_edited |>
                                    dplyr::mutate(table = "EBS shelf-extended")) |>
                 dplyr::mutate(table = factor(table,
                                              levels = table_order)),
               mapping = aes(x = min_depth_fm,
                             xend = max_depth_fm,
                             y = scope_fm,
                             yend = scope_fm,
                             color = table,
                             linetype = table),
               linewidth = 1
               ) +
  scale_color_manual(values = c("EBS slope" = "orange",
                                "EBS shelf" = "#CC79A7",
                                "GOA/AI" = "#009E73",
                                "EBS slope-modified" = "#000000", 
                                "EBS shelf-extended" = "#56B4E9")) +
  scale_x_continuous(name = "Depth (fm)") +
  scale_y_continuous(name = "Scope (fm)") +
  ggtitle("Modified Scope Tables 2024") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.2))


plot_slope_table_m <- ggplot() +
  geom_segment(data = 
                 dplyr::bind_rows(modified_slope_edited |>
                                    dplyr::mutate(table = "EBS slope-modified"),
                                  modified_shelf_edited |>
                                    dplyr::mutate(table = "EBS shelf-extended")) |>
                 dplyr::mutate(table = factor(table,
                                              levels = table_order)),
               mapping = aes(x = min_depth_m,
                             xend = max_depth_m,
                             y = scope_m,
                             yend = scope_m,
                             color = table,
                             linetype = table),
               linewidth = 1
  ) +
  scale_color_manual(values = c("EBS slope" = "orange",
                                "EBS shelf" = "#CC79A7",
                                "GOA/AI" = "#009E73",
                                "EBS slope-modified" = "#000000", 
                                "EBS shelf-extended" = "#56B4E9")) +
  scale_x_continuous(name = "Depth (m)") +
  scale_y_continuous(name = "Scope (m)") +
  ggtitle("Modified Scope Tables 2024") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.2))

