# Tables and figures comparing fits to Somerton data


load(here::here(output_dir, "results_2002.rda"))


ols_perf <- 
  lapply(
    results_2002,
    FUN = function(x) { x$loocv_table}) |>
  do.call(what = rbind) |>
  dplyr::filter(model_name == "ols1") |>
  dplyr::select(common_name, method, pbias, rmse) |>
  dplyr::mutate(
    rmse = trimws(format(round(rmse), big.mark = ",")),
    pbias = trimws(format(round(pbias, 1), nsmall = 1))
    )


ols_par <-
  lapply(
    results_2002,
    FUN = function(x) {
      out <- x$ols_results$fpc
      out$common_name <- x$common_name
      out
    }
  ) |>
  do.call(what = rbind) |>
  dplyr::mutate(
    log_ratio = format(round(log_ratio, 3), nsmall = 3),
    ratio = format(round(ratio, 3), nsmall = 3),
    ratio_lci = format(round(ratio_lci, 3), nsmall = 3),
    ratio_uci = format(round(ratio_uci, 3), nsmall = 3),
    var = format(round(var, 3), nsmall = 3),
    ratio = paste0(ratio, " (", ratio_lci, "-", ratio_uci, ")")
  )

ols_table <- 
  dplyr::inner_join(
    ols_par,
    ols_perf
  )

dplyr::select(ols_table, common_name, method, ratio, rmse, pbias) |>
  as.data.frame() |>
  xlsx::write.xlsx(
    file = here::here("analysis", "somerton_2002", "plots", "ols_comparison_table.xlsx"),
    row.names = FALSE
  )


oos_results_2002 <- 
  lapply(
  results_2002,
  FUN = function(x) { x$oos_table}) |>
  do.call(what = rbind) |>
  dplyr::filter(model_name == "ols1")


p_oos_2002 <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    data = oos_results_2002,
    mapping = aes(x = common_name, y = pbias, color = method),
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = oos_results_2002,
    mapping = aes(x = common_name, ymin = pbias_lci, ymax = pbias_uci, color = method),
    position = position_dodge(width = 0.5),
    width = 0
  ) +
  scale_y_continuous(name = "PBIAS (%)") +
  scale_color_manual(values = c("grey", "black")) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.22,0.88),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8),    
    axis.title.y = element_text(size = 8.5),
    axis.title.x = element_blank()
  )


png(filename = here::here("analysis", "somerton_2002", "plots", "pbias_2002_vs_other.png"),
    width = 80,
    height = 80,
    units = "mm",
    res = 600)
print(p_oos_2002)
dev.off()
