library(trawlmetrics)
library(cowplot)

ebs_nbs_catch <- trawlmetrics::bts_geom |>
  dplyr::filter(SURVEY_DEFINITION_ID %in% c(98, 143))

p1 <- 
  ggplot() +
  geom_histogram(data = ebs_nbs_catch,
                 mapping = aes(x = TOTAL_WEIGHT_KG)) +
  geom_vline(xintercept = 900, linetype = 2) +
  scale_x_log10(name = "Total weight (kg)", limits = c(10, 17000)) +
  scale_y_continuous(name = "Hauls (#)", limits = c(0, 2000), expand = c(0, 0)) +
  facet_wrap(~"EBS/NBS survey catch (30 min. tows)") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18))

p2 <- 
  ggplot() +
  geom_histogram(data = ebs_nbs_catch,
                 mapping = aes(x = TOTAL_WEIGHT_KG/2)) +
  geom_vline(xintercept = 900, linetype = 2) +
  scale_x_log10(name = "Total weight (kg)", limits = c(10, 17000)) +
  scale_y_continuous(name = "Hauls (#)", limits = c(0, 2000), expand = c(0, 0)) +
  facet_wrap(~"Projected (15 min. tows)") +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18))

p_grid <- 
  cowplot::plot_grid(
    p1, p2, align = "hv", nrow = 2
  )

ragg::agg_png(
  filename = here::here("analysis", "15_30", "plots", "historical_total_catch.png"),
  width = 8, 
  height = 6,
  units = "in",
  res = 300
)
print(p_grid)
dev.off()
