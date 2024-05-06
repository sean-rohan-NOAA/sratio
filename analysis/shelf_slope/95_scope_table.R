library(sratio)

# 2023 Scope Values
scope_table_2023 <- readRDS(here::here("analysis", "shelf_slope", "data", "ss_haul.rds")) |> 
  dplyr::select(dplyr::all_of(c("BOTTOM_DEPTH", "GEAR", "WIRE_LENGTH"))) |>
  dplyr::inner_join(data.frame(GEAR = c(44, 172), GEAR_NAME = c("83-112", "PNE")))


# EBS shelf and slope standard scope tables
scope_table_standard <- readxl::read_xlsx(
  here::here("analysis", "shelf_slope", "docs", "shelf_slope_table.xlsx")
  )


# Plot 2023 scope table
ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "2023_depth_vs_scope.png"),
              width = 6, height = 4, units = "in", res = 180)
print(
  ggplot() +
    geom_point(data = scope_table_2023,
               mapping = aes(x = BOTTOM_DEPTH, 
                             y = WIRE_LENGTH, 
                             color = GEAR_NAME, 
                             shape = GEAR_NAME)) +
    scale_x_continuous(name = "Bottom Depth (m)") +
    scale_y_continuous(name = "Wire Length (m)") +
    ggtitle(label = "2023 Shelf/Slope Scope") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
)
dev.off()




