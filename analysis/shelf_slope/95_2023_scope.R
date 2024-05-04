library(sratio)

scope_table <- readRDS(here::here("analysis", "shelf_slope", "data", "ss_haul.rds")) |> 
  dplyr::select(dplyr::all_of(c("BOTTOM_DEPTH", "GEAR", "WIRE_LENGTH"))) |>
  dplyr::inner_join(data.frame(GEAR = c(44, 172), GEAR_NAME = c("83-112", "PNE")))


ragg::agg_png(filename = here::here("analysis", "shelf_slope", "plots", "2023_depth_vs_scope.png"),
              width = 6, height = 4, units = "in", res = 180)
print(
  ggplot() +
    geom_point(data = scope_table,
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
