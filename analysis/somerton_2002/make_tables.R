library(xlsx)
library(ggplot2)
library(ggthemes)
library(glmmTMB)
library(scales)
library(dplyr)
library(ggpp)


project_dir <- here::here("analysis", "somerton_2002")

output_dir <- here::here("analysis", "somerton_2002", "output")
                         
# Fits to 1998 data with both sexes in the same model (Somerton et al. (2002); validation with 1995+2021-2024
load(here::here(output_dir, "results_2002.rda"))

# Models fitted to 1998; validation with 1995+2021-2024
load(here::here(output_dir, "results_95_2124.rda"))

# Models fitted to all data; no two-fold cross validation
load(here::here(output_dir, "results_all.rda"))

# Models fitted to 1998 data, validation with 1995+2021-2024
load(here::here(output_dir, "results_1998.rda"))


# Make RMSE tables by combining LOOCV, OOS (two-fold CV), and diagnostic (AIC) tables.
tables_rmse <- 
  function(results_obj, save_dir) {
  rmse_output <- 
    lapply(
      results_obj,
      FUN = function(x) {
        
        subset_name <- x$subset_name
        common_name <- x$common_name
        
        fname <- 
          here::here(save_dir, "plots", paste0(subset_name, "_fits"),  paste0("rmse_table_", subset_name, "_", common_name, ".xlsx"))
        
        output <-
          dplyr::inner_join(
            x[['loocv_table']],
            x[['aic_table']]
          )|>
          dplyr::select(
            common_name,
            method,
            rmse,
            pbias,
            convergence,
            pdhess,
            max_gradient, 
            pass_check
          ) |>
          dplyr::arrange(
            rmse
          ) 
        
        if(is(x[['oos_table']], "data.frame")) {
          
          oos <- x[['oos_table']]
          
          names(oos)[grepl(pattern = "rmse", x = names(oos))] <- 
            paste0(
              "oos_", names(oos)[grepl(pattern = "rmse", x = names(oos))]
            )
          
          names(oos)[grepl(pattern = "pbias", x = names(oos))] <- 
            paste0(
              "oos_", names(oos)[grepl(pattern = "pbias", x = names(oos))]
            )
          
          oos <- oos |>
            dplyr::select(-best)
          
          oos$contrast_name <- x$contrast_name
          
          output <- 
            dplyr::right_join(
            output, oos
          )
          
        }
        
        write.xlsx2(
          x = output, 
          file = fname, 
          row.names = FALSE,
          sheetName = common_name
        )
        
        output
      }
      
    )
  
  return(rmse_output)
  }

combine_fits <- function(results_obj) {
  rmse_output <- 
    lapply(
      results_obj,
      FUN = function(x) {
        
        subset_name <- x$subset_name
        common_name <- x$common_name
        
        output <- x[['fit_table']]
          dplyr::inner_join(
            
          )
      }
    )
  
}
  

# Make final RMSE table ----

response_type <- 
  data.frame(
    method = c("OLS median", "OLS mean", paste0("lognormal", 1:5),  paste0("ccr_beta", 1:15), paste0("ccr_bin", 1:5), paste0("bin", 1:5), paste0("bb", 1:15), paste0("pois", 1:8), paste0("nb", 1:24)),
    type = c(rep("Ratio", 7), rep("CCR", 20), rep("Proportion", 20), rep("Count", 32)),
    type_abbv = c(rep("Ratio", 7), rep("CCR", 20), rep("Prop", 20), rep("Count", 32))
  )

rmse_all <-  
  tables_rmse(
    results_obj = results_all,
    save_dir = project_dir 
  )

rmse_1998 <-  
  tables_rmse(
    results_obj = results_1998,
    save_dir = project_dir 
  )

rmse_95_2124 <-  
  tables_rmse(
    results_obj = results_95_2124,
    save_dir = project_dir
  )


# RMSE table for converged models -- model with the lowest RMSE for each type

best_rmse <- 
  do.call(what = rbind, rmse_all) |>
  dplyr::inner_join(response_type ) |>
  dplyr::filter(pass_check == TRUE) |>
  dplyr::select(common_name, type, method, rmse, pbias) |>
  dplyr::group_by(common_name, type) |>
  dplyr::filter(rmse == min(rmse)) |>
  dplyr::ungroup() |>
  dplyr::group_by(common_name) |>
  dplyr::mutate(w = sqrt(rmse)^-2/ sum(sqrt(rmse)^-2))


write.csv(
  x = dplyr::arrange(best_rmse, common_name, rmse) |>
    dplyr::mutate(rmse = round(rmse),
                  pbias = format(round(pbias, 1), nsmall = 1),
                  w = format(round(w, 3), nsmall = 1)),
  file = here::here(project_dir, "plots", "lowest_rmse_models_table.csv"),
  row.names = FALSE
)


do.call(what = rbind, rmse_all) |>
  dplyr::filter(method == "OLS mean")

# 2CV prediction performance for the 'best' models in each category

oos_performance <- 
  dplyr::bind_rows(
    best_rmse |>
      dplyr::inner_join(
        do.call(what = rbind, rmse_95_2124) |>
          dplyr::inner_join(response_type) |>
          dplyr::select(common_name, type, type_abbv, method, oos_rmse, oos_rmse_lci, oos_rmse_uci, oos_pbias, oos_pbias_lci, oos_pbias_uci, contrast_name)
      ),
    best_rmse |>
      dplyr::inner_join(
        do.call(what = rbind, rmse_1998) |>
          dplyr::inner_join(response_type) |>
          dplyr::select(common_name, type, type_abbv, method, oos_rmse, oos_rmse_lci, oos_rmse_uci, oos_pbias, oos_pbias_lci, oos_pbias_uci, contrast_name)
      )
  ) |>
  dplyr::arrange(
    common_name, rmse
  ) |>
  dplyr::mutate(
    validation_subset = ifelse(contrast_name == "1998FitVsOther", "1995&2021-2024", "1998")
  )

p_oos_pbias <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(
    data = oos_performance,
    mapping = aes(x = type_abbv, y = oos_pbias, color = validation_subset),
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = oos_performance,
    mapping = aes(x = type_abbv, ymin = oos_pbias_lci, ymax = oos_pbias_uci, color = validation_subset),
    width = 0,
    position = position_dodge(width = 0.5)
  ) +
  scale_y_continuous(name = "PBIAS (%)",
                     limits = c(-100, 150),
                     oob = oob_squish,
                     expand = c(0,0)
  ) +
  scale_x_discrete(name = "Response") +
  scale_color_manual(name = "Validation dataset", values = c("#40B0A6", "#5D3A9B")) +
  facet_wrap(~common_name, ncol = 1) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8.5),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.5, "mm")
  )

# Compare models that converged

  

compare_perf_among_models <- 
  do.call(what = rbind, rmse_all) |>
  dplyr::filter(pass_check == TRUE) |>
  dplyr::group_by(
    common_name
  ) |>
  dplyr::mutate(
    rel_rmse = rmse/min(rmse)
  ) |>
  dplyr::inner_join(
    response_type
  ) |>
  dplyr::mutate(
    method = factor(method, levels = response_type$method)
  )
# ggplot() +
#   geom_point(data = compare_perf_among_models,
#              mapping = aes(y = method, x = rel_rmse, color = common_name),
#              shape = 21) +
#   scale_color_viridis_d(name = "Common name", option = "H") +
#   scale_x_continuous(name = "RMSE/min(RMSE)") +
#   theme_bw() +
#   theme(legend.position = "inside",
#         legend.position.inside = c(0.8, 0.8),
#         axis.title.y = element_blank())

p_bias_method <- 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = c(-10,10), linetype = 3) +
  geom_point(data = compare_perf_among_models,
             mapping = aes(y = method, x = pbias, color = type, shape = common_name)) +
  geom_text(
    data = compare_perf_among_models |>
      dplyr::group_by(method, type) |>
      dplyr::summarise(n = n()),
    mapping = aes(x = -54, y = method, label = n, color = type),
    size = 3
  ) +
  scale_color_tableau(name = "Response type") +
  scale_x_continuous(name = "PBIAS (%)", breaks = seq(-50, 100, 25)) +
  scale_shape(name = "Species/sex", solid = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 9),
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(3.5, "mm"))

png(filename = here::here("analysis", "somerton_2002", "plots", "PBIAS_by_method_type_species.png"),
    width = 169,
    height = 169,
    units = "mm",
    res = 600)
print(p_bias_method)
dev.off()


# ggplot() +
#   geom_text(data = compare_perf_among_models,
#              mapping = aes(x = pbias, y = rel_rmse, color = type, label = method)) +
#   scale_x_continuous(name = "|PBIAS| (%)") +
#   scale_y_continuous(name = "RMSE") +
#   theme_bw() +
#   theme(legend.position = "inside",
#         legend.position.inside = c(0.8, 0.8),
#         axis.title.y = element_blank())


# Supplementary table of all RMSEs for models that converged ----
all_rmse <- 
  do.call(what = rbind, rmse_all) |>
  dplyr::inner_join(response_type ) |>
  dplyr::filter(pass_check == TRUE) |>
  dplyr::select(common_name, type, method, rmse, pbias) |>
  dplyr::left_join(
    do.call(what = rbind, rmse_95_2124) |>
      dplyr::inner_join(response_type ) |>
      dplyr::filter(pass_check == TRUE) |>
      dplyr::select(common_name, type, method, rmse_1998 = rmse, pbias_1998 = pbias)
  ) |>
  dplyr::left_join(
    do.call(what = rbind, rmse_1998) |>
      dplyr::inner_join(response_type) |>
      dplyr::filter(pass_check == TRUE) |>
      dplyr::select(common_name, type, method, rmse_95_2124 = rmse, pbias_95_2124 = pbias)
  ) |>
  dplyr::arrange(
    common_name, rmse
  )


# Plot fits for the best-fit model in each category

best_fits <- 
  lapply(results_all, 
         FUN = function(results) {
           results$fit_table
         }
           ) |>
  do.call(what = rbind) |>
  dplyr::mutate(method = ifelse(model_name == "ols1", method, model_name)) |>
  dplyr::inner_join(
    best_rmse
  ) |>
  dplyr::left_join(
    response_type
  )

cpue_dat <-
  lapply(results_all, 
         FUN = function(results) {
           results$dat
         }
  ) |>
  do.call(what = rbind) |>
  dplyr::mutate(
    common_name = paste0(common_name, " (", ifelse(SEX == 'M', "male", "female"), ")")
  )




p_obs_fit <- 
  ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(
    data = cpue_dat,
    mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_30),
    color = "grey50",
    size = 1,
    shape = 21
  ) +
  geom_ribbon(
    data = best_fits, 
    mapping = aes(x = CPUE_NO_KM2_15, ymin = fit_lwr, ymax = fit_upr, fill = type),
    alpha = 0.1
  ) +
  geom_path(data = best_fits, 
            mapping = aes(x = CPUE_NO_KM2_15, y = fit, color = type)) +
  ggpp::geom_text_npc(
    data = somerton_estimates,
    mapping = aes(npcx = 0.02, npcy = .95, label = common_name),
    size = 2.7) +
  facet_wrap(~common_name, scales = "free", ncol = 1) +
  scale_x_log10(name = expression(CPUE[15]*' (#/'*km^2*')')) +
  scale_y_log10(name = expression(CPUE[30]*' (#/'*km^2*')')) +
  scale_color_colorblind(name = "Type") +
  scale_fill_colorblind(name = "Type") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8),    
    axis.title = element_text(size = 8.5),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.5, "mm")
  )+ # Set legend to horizontal position
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


# Somerton FPC fits
somerton_estimates <-
  data.frame(
    common_name = c("snow crab (male)", "snow crab (female)", "red king crab (male)", "red king crab (female)", "Tanner crab (male)", "Tanner crab (female)"),
    fpc = c(1.78, 1.78, 1.244, 1.244, 1.68, 1.68)
  )

p_fishing_power <- 
  ggplot() +
  geom_abline(slope = 0, intercept = 1, linetype = 2) +
  geom_abline(
    data = somerton_estimates,
    mapping = aes(slope = 0, intercept = fpc),
    linetype = 3, color = "red"
  ) +
  geom_rug(
    data = cpue_dat,
    mapping = aes(x = CPUE_NO_KM2_15),
    color = "grey50",
    length = unit(1, "mm")
  ) +
  geom_ribbon(
    data = best_fits,
    mapping = aes(x = CPUE_NO_KM2_15, ymin = CPUE_NO_KM2_15/fit_lwr, ymax = CPUE_NO_KM2_15/fit_upr, fill = type),
    alpha = 0.2
  ) +
  geom_path(data = dplyr::arrange(best_fits, CPUE_NO_KM2_15), 
            mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_15/fit, color = type)) +
  facet_wrap(~common_name, scales = "free_x", ncol = 1) +
  scale_x_log10(name = expression(CPUE[15]*' (#/'*km^2*')')) +
  scale_y_continuous(name = expression('Fishing power ratio ('*q[15]/q[30]*')'), expand = c(0,0), oob = scales::oob_censor) +
  scale_color_colorblind(name = "Type") +
  scale_fill_colorblind(name = "Type") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8),    
    axis.title = element_text(size = 8.5),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.5, "mm")
  )

p_multipanel_fpc <- 
  cowplot::plot_grid(
    p_obs_fit,
    p_fishing_power,
    p_oos_pbias,
    align = "hv",
    ncol = 3
  )

png(
  filename = here::here("analysis", "somerton_2002", "plots", "fpc_fit_all_years.png"), 
  width = 169, 
  height = 169,
  units = "mm",
  res = 600
)
print(p_multipanel_fpc)
dev.off()


## 

p_fishing_power_multipanel <-
  ggplot() +
  geom_abline(slope = 0, intercept = 1, linetype = 2) +
  geom_abline(
    data = somerton_estimates,
    mapping = aes(slope = 0, intercept = fpc),
    linetype = 3, color = "red"
  ) +
  geom_ribbon(
    data = best_fits,
    mapping = aes(x = CPUE_NO_KM2_15, ymin = CPUE_NO_KM2_15/fit_lwr, ymax = CPUE_NO_KM2_15/fit_upr, fill = type),
    alpha = 0.2
  ) +
  geom_path(data = dplyr::arrange(best_fits, CPUE_NO_KM2_15), 
            mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_15/fit, color = type)) +
  geom_rug(
    data = cpue_dat,
    mapping = aes(x = CPUE_NO_KM2_15),
    color = "grey50",
    length = unit(1, "mm")
  ) +
  facet_grid(type~common_name, scales = "free_x") +
  scale_x_log10(name = expression(CPUE[15]*' (#/'*km^2*')')) +
  scale_y_continuous(name = expression('Fishing power ratio ('*q[15]/q[30]*')'), expand = c(0,0), oob = scales::oob_keep) +
  scale_color_colorblind(name = "Type") +
  scale_fill_colorblind(name = "Type") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text.x = element_text(size = 5.75),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),    
    axis.text.y = element_text(size = 8),    
    axis.title = element_text(size = 8.5),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.5, "mm")
  )

png(
  filename = here::here("analysis", "somerton_2002", "plots", "fpc_fit_all_years_multipanel.png"), 
  width = 169, 
  height = 169,
  units = "mm",
  res = 600
)
print(p_fishing_power_multipanel)
dev.off()
