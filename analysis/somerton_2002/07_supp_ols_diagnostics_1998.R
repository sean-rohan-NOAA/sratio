# Make OLS residual diagnostic plots for models fitted to 2002 data

library(glmmTMB)
library(nortest)
library(ggplot2)
library(dplyr)
library(ggthemes)

load(here::here("analysis", "somerton_2002", "output", "results_2002.rda"))

unique_name <- names(results_2002)

for(uu in unique_name) {
  
  sel_res <- results_2002[[uu]]
  
  model_data <- fortify(sel_res$ols_results$best_model)
  
  p_hist <- 
    ggplot(data = model_data) +
    geom_histogram(
      mapping = aes(x = .stdresid),
      bins = 15
    ) +
    scale_x_continuous(name = "Std. residuals") +
    scale_y_continuous(name = "Frequency") +
    theme_bw()
  
  kurt_val <- format(sel_res$ols_results$kurtosis, digits = 2, nsmall = 2)
  ad_stat  <- format(unname(sel_res$ols_results$anderson_darling$statistic), digits = 1, nsmall = 1)
  p_val    <- sel_res$ols_results$anderson_darling$p.value
  p_text   <- if(p_val < 0.001) "< 0.001" else paste("==", format(p_val, digits = 1, nsmall = 3))
  
  qq_lab <- paste0(
    "'Kurtosis: '*G[2] == ", kurt_val, 
    " ~~~~~ 'A-D: p ' ", p_text
  )
  
  p_qq <-
    ggplot(
      data = model_data,
      aes(sample = .stdresid)
    ) +
    stat_qq() +
    stat_qq_line(linetype = 2) +
    annotate(geom = "text", x = -Inf, y = Inf, label = qq_lab, hjust = -0.05, vjust = 1.5, parse = TRUE) +
    labs(x = "Theoretical Quantiles",
         y = "Std. Residuals") +
    theme_bw()
  
  
  
  png(
    filename = here::here("analysis", "somerton_2002", "plots", "2002_fits", paste0("OLS_diagnostics_", uu, ".png")),
    width = 169,
    height = 169,
    units = "mm",
    res = 600
  )
  print(cowplot::plot_grid(
    sel_res$p_heteroskedasticity, 
    cowplot::plot_grid(p_hist, p_qq),
    nrow = 2
  ))
  dev.off()
  
}



