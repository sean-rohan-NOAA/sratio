library(xlsx)
library(ggplot2)
library(ggthemes)
library(glmmTMB)
library(scales)
library(dplyr)
library(ggpp)
library(stringr)


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


# Setup CPUE data

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

determine_poly_order <- function(formula_input) {
  
  # Helper function to parse a single formula (character or formula object)
  parse_single <- function(f) {
    # 1. Convert formula objects to a character string safely
    if (inherits(f, c("formula", "Formula"))) {
      f_str <- paste(deparse(f), collapse = " ")
    } else {
      f_str <- as.character(f)
    }
    
    # 2. Clean up double-tilde formats (e.g., stripping out the trailing " ~1")
    f_str <- trimws(f_str)
    parts <- strsplit(f_str, "\\s*~\\s*")[[1]]
    if (length(parts) > 2) {
      if (parts[1] == "") {
        f_str <- paste0("~ ", parts[2]) # One-sided formula
      } else {
        f_str <- paste0(parts[1], " ~ ", parts[2]) # Two-sided formula
      }
    }
    
    # 3. Convert back to R formula and parse terms
    f_obj <- as.formula(f_str)
    t_obj <- terms(f_obj)
    
    # 4. Check for the presence of an intercept
    has_intercept <- attr(t_obj, "intercept") == 1
    
    # 5. Extract term labels and find the maximum polynomial degree
    term_labels <- attr(t_obj, "term.labels")
    
    if (length(term_labels) == 0) {
      max_deg <- 0 # Intercept-only model
    } else {
      degrees <- sapply(term_labels, function(label) {
        # Match poly(x, degree)
        poly_match <- regmatches(label, regexec("poly\\(.*,\\s*([0-9]+)", label))[[1]]
        if (length(poly_match) > 1) {
          return(as.numeric(poly_match[2]))
        }
        
        # Match I(x^degree) or x^degree
        pow_match <- regmatches(label, regexec("\\^\\s*([0-9]+)", label))[[1]]
        if (length(pow_match) > 1) {
          return(as.numeric(pow_match[2]))
        }
        
        # Standard linear terms default to degree 1
        return(1)
      })
      max_deg <- max(degrees)
    }
    
    # 6. Format output based on intercept and degree
    if (has_intercept) {
      return(as.character(max_deg))
    } else {
      return(paste0("0+", max_deg))
    }
  }
  
  # Vectorize to support character vectors or lists of formulas
  if (inherits(formula_input, c("formula", "Formula"))) {
    return(parse_single(formula_input))
  } else if (is.list(formula_input) || length(formula_input) > 1) {
    return(sapply(formula_input, parse_single, USE.NAMES = FALSE))
  } else {
    return(parse_single(formula_input))
  }
}

# Model name and polynomial order for mean and variance formulas. Zero-intercept models are noted with a leading '0+', i.e., 0+1 for linear 


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
        
        aic_table <- x[['aic_table']]
        aic_table$disp[is.na(aic_table$disp)] <- "~1" # Fill OLS order
        
        # Get polynomial order for each model
        aic_table$order_mean <- determine_poly_order(aic_table$formula)
        aic_table$order_disp <- determine_poly_order(aic_table$disp)

        
        output <-
          dplyr::inner_join(
            x[['loocv_table']],
            aic_table
          ) |>
          dplyr::select(
            common_name,
            method,
            rmse,
            pbias,
            convergence,
            pdhess,
            max_gradient, 
            pass_check,
            order_mean,
            order_disp
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
  


# Make RMSE tables ----

response_type <- 
  data.frame(
    method = c("OLS median", "OLS mean", paste0("LN", 1:5),  paste0("CCR_BB", 1:15), paste0("CCR_BI", 1:5), paste0("BI", 1:5), paste0("BB", 1:15), paste0("PO", 1:8), paste0("NB", 1:24)),
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


# Review fits --------------------------------------------------------------------------------------

# Compare models that converged
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
    method = factor(method, levels = response_type$method),
    method_order = paste0(method, "(", order_mean, "; ", order_disp,")")
  )

used_models <- left_join(
  response_type,
  compare_perf_among_models |>
    dplyr::ungroup() |>
    dplyr::select(method, method_order) |>
    unique()) |>
  dplyr::filter(!is.na(method_order))

compare_perf_among_models$method_order <- factor(compare_perf_among_models$method_order, levels = used_models$method_order)
  
xlsx::write.xlsx(all_rmse, file = here::here("analysis", "somerton_2002", "plots", "rmse_all_converged.xlsx"), row.names = FALSE)

p_bias_method <- 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = c(-10,10), linetype = 3) +
  geom_point(data = compare_perf_among_models,
             mapping = aes(y = method_order, x = pbias, color = type, shape = common_name)) +
  geom_text(
    data = compare_perf_among_models |>
      dplyr::group_by(method_order, type) |>
      dplyr::summarise(n = n()),
    mapping = aes(x = -54, y = method_order, label = n, color = type),
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
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.3),
        legend.key.height = unit(3.5, "mm"))

png(filename = here::here("analysis", "somerton_2002", "plots", "PBIAS_by_method_type_species.png"),
    width = 169,
    height = 169,
    units = "mm",
    res = 600)
print(p_bias_method)
dev.off()

do.call(what = rbind, rmse_all) |>
  dplyr::filter(method == "OLS mean")



# Somerton FPC fits
somerton_estimates <-
  data.frame(
    common_name = c("snow crab (male)", "snow crab (female)", "red king crab (male)", "red king crab (female)", "Tanner crab (male)", "Tanner crab (female)"),
    fpc = c(1.78, 1.78, 1.244, 1.244, 1.68, 1.68)
  )

# Plot converged models for each species 

## Show fits for each converged model
all_fits <- 
  lapply(
    X = names(results_all),
    FUN = function(spname) {
      
      results_all[[spname]]$fit_table |>
        dplyr::mutate(
          method = ifelse(is.na(method), model_name, method)) |>
        dplyr::inner_join(all_rmse) |>
        dplyr::inner_join(response_type)
    }
  )

species_model_labs <- 
  all_rmse |>
  dplyr::mutate(
    rmse = round(rmse),
    pbias = format(round(pbias, 1), nsmall = 1),
    label = paste0("RMSE=",rmse, "\nPBIAS=",trimws(pbias))
  )

p_all_models_by_species <- vector(mode = "list", length = length(all_fits))

for(hh in 1:length(all_fits)) {
  
  sel_common_name <- all_fits[[hh]]$common_name[1]
  
  p_all_models_by_species[[hh]] <- 
    ggplot() +
    geom_abline(slope =1 , intercept = 0, linetype = 2) +
    geom_point(
      data = dplyr::filter(cpue_dat, common_name == sel_common_name),
      mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_30),
      color = "grey50",
      size = 0.7
    ) +
    geom_ribbon(
      data = all_fits[[hh]],
      mapping = aes(x = CPUE_NO_KM2_15, ymin = fit_lwr, ymax = fit_upr, fill = type),
      alpha = 0.3
    ) +
    geom_path(
      data = all_fits[[hh]],
      mapping = aes(x = CPUE_NO_KM2_15, y = fit, color = type)
    ) +
    ggpp::geom_text_npc(
      data = dplyr::filter(species_model_labs, common_name == sel_common_name),
      mapping = aes(npcx = 0.02, npcy = .97, label = method),
      size = 2) +
    ggpp::geom_text_npc(
      data = dplyr::filter(species_model_labs, common_name == sel_common_name),
      mapping = aes(npcx = 0.97, npcy = .05, label = label ),
      size = 2) +
    facet_wrap(~factor(method, levels = response_type$method)) +
    scale_x_log10(name = expression(CPUE[15]*' (#/'*km^2*')')) +
    scale_y_log10(name = expression(CPUE[30]*' (#/'*km^2*')')) +
    scale_color_colorblind(name = "Type") +
    scale_fill_colorblind(name = "Type") +
    ggtitle(paste(sel_common_name, "fits")) +
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
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
  
  png(filename = here::here("analysis", "somerton_2002", "plots", paste0("fits_obs_by_species_method_", sel_common_name, ".png")),
      width = 254,
      height = 150,
      units = "mm",
      res = 600)
  print(p_all_models_by_species[[hh]])
  dev.off()
  
}


# Select best models after diagnostic checks (outside of R) ----------------------------------------

# Models were screened for goodness-of-fit based on DHARMa residuals (sometimes QQ plots when 
# DHARMa residuals failed to generate (ccr_bin models). Spreadsheet includes an explanation of
# why each model was ultimately selected.

final_models <- 
  xlsx::read.xlsx(
  file = here::here("analysis", "somerton_2002", "data", "best_models.xlsx"),
  sheetIndex = 1
) |>
  dplyr::inner_join(
    dplyr::select(all_rmse, method, common_name, rmse, pbias)
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(common_name) |>
  dplyr::mutate(w = sqrt(rmse)^-2/ sum(sqrt(rmse)^-2)) |>
  dplyr::arrange(
    desc(common_name), type
  )

xlsx::write.xlsx(
  final_models |>
    dplyr::mutate(
      rmse = round(rmse),
      pbias = format(round(pbias, 1), nsmall = 1),
      w = format(round(w, 3), nsmall = 1)
    ) |>
    as.data.frame(),
  file = here::here("analysis", "somerton_2002", "plots", "final_models.xlsx"),
  row.names = FALSE
)


# 2CV prediction performance for the 'best' models in each category
oos_performance <- 
  dplyr::bind_rows(
    final_models |>
      dplyr::inner_join(
        do.call(what = rbind, rmse_95_2124) |>
          dplyr::inner_join(response_type) |>
          dplyr::select(common_name, type, type_abbv, method, oos_rmse, oos_rmse_lci, oos_rmse_uci, oos_pbias, oos_pbias_lci, oos_pbias_uci, contrast_name)
      ),
    final_models |>
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
    final_models
  ) |>
  dplyr::left_join(
    response_type
  )


best_fits |>
  dplyr::group_by(common_name, type, model_name) |>
  dplyr::summarise(fit_mean = round(mean(CPUE_NO_KM2_15/fit), 2),
                   fit_min = round(min(CPUE_NO_KM2_15/fit), 2),
                   fit_max = round(max(CPUE_NO_KM2_15/fit), 2),
                   fit_lwr_mean = round(mean(CPUE_NO_KM2_15/fit_lwr), 2),
                   fit_lwr_min = round(min(CPUE_NO_KM2_15/fit_lwr), 2), 
                   fit_lwr_max = round(max(CPUE_NO_KM2_15/fit_lwr), 2), 
                   fit_upr_mean = round(mean(CPUE_NO_KM2_15/fit_upr), 2),
                   fit_upr_min = round(min(CPUE_NO_KM2_15/fit_upr), 2), 
                   fit_upr_max = round(max(CPUE_NO_KM2_15/fit_upr), 2)
  ) |>
  as.data.frame() |>
  xlsx::write.xlsx(
    file = here::here("analysis", "somerton_2002", "plots", "mean_fpc_by_species.xlsx"),
    row.names = FALSE
  )


# How common were density effects among best-fit models?
model_formulas <- 
  lapply(results_all, 
         FUN = function(results) {
           dplyr::select(results$aic_table, common_name, model_name, formula, disp)
         }
  ) |>
  do.call(what = rbind) |>
  dplyr::right_join(
    dplyr::mutate(final_models, model_name = ifelse(method %in% c("OLS median", "OLS mean"), "ols1", method))
  )

# True = no effect on the mean
table(grepl(pattern = " ~ 1", x = model_formulas$formula))
table(grepl(pattern = " ~ 1", x = model_formulas$formula), model_formulas$common_name)

# True = no effect on thevariance
table(grepl(pattern = "~1", x = model_formulas$disp[model_formulas$type != "Ratio"]), model_formulas$type[model_formulas$type != "Ratio"])
table(grepl(pattern = "~1", x = model_formulas$disp[model_formulas$type != "Ratio"]), 
      model_formulas$common_name[model_formulas$type != "Ratio"])

dplyr::filter(model_formulas, common_name == "red king crab (male)")

# Observed versus fitted
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
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))


rmse_labs <- 
  final_models |>
  dplyr::mutate(
    rmse = round(rmse),
    pbias = format(round(pbias, 1), nsmall = 1),
    w = format(round(w, 3), nsmall = 1),
    label = paste0(method, "\nRMSE=",rmse, "\nPBIAS=",trimws(pbias), "%")
  )

p_obs_fit_multipanel <- 
  ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(
    data = cpue_dat,
    mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_30, color = factor(YEAR), shape = factor(YEAR)),
    size = 1
  ) +
  geom_ribbon(
    data = best_fits, 
    mapping = aes(x = CPUE_NO_KM2_15, ymin = fit_lwr, ymax = fit_upr,),
    alpha = 0.1
  ) +
  geom_path(data = best_fits, 
            mapping = aes(x = CPUE_NO_KM2_15, y = fit)) +
  geom_path(data = dplyr::inner_join(best_fits, somerton_estimates),
            mapping = aes(x = CPUE_NO_KM2_15, y = CPUE_NO_KM2_15/fpc),
            color = "red",
            linetype = 3) +
  ggpp::geom_text_npc(
    data = rmse_labs,
    mapping = aes(npcx = 0.02, npcy = .97, label = label),
    size = 2.5) +
  facet_grid(type~common_name, scales = "free_x") +
  scale_x_log10(name = expression(CPUE[15]*' (#/'*km^2*')')) +
  scale_y_log10(name = expression(CPUE[30]*' (#/'*km^2*')')) +
  scale_color_colorblind(name = "Year") +
  scale_shape(name = "Year") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    # legend.title.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text = element_text(size = 8),    
    axis.title = element_text(size = 8.5),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.5, "mm")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

png(filename = here::here("analysis", "somerton_2002", "plots", "fits_obs_by_year_method.png"),
    width = 254,
    height = 150,
    units = "mm",
    res = 600)
print(p_obs_fit_multipanel)
dev.off()


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
    # p_oos_pbias,
    align = "hv",
    ncol = 2
    # ncol = 3
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
  ggpp::geom_text_npc(
    data = rmse_labs,
    mapping = aes(npcx = 0.02, npcy = .97, label = label),
    size = 2.5) +
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
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text.x = element_text(size = 6.75),
    strip.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),    
    axis.text.y = element_text(size = 8),    
    axis.title = element_text(size = 8.5),
    legend.key.height = unit(2.5, "mm"),
    legend.key.width = unit(2.5, "mm")
  )

png(
  filename = here::here("analysis", "somerton_2002", "plots", "fpc_fit_all_years_multipanel.png"), 
  width = 254, 
  height = 150,
  units = "mm",
  res = 600
)
print(p_fishing_power_multipanel)
dev.off()
