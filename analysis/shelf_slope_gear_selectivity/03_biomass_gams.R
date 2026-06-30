# library(sratio)
library(glmmTMB)
library(DHARMa)

catch_data <- readRDS(
  here::here("analysis", "shelf_slope_gear_selectivity", "data", "catch_data.rds")
)

catch_data_wide <- readRDS(
  here::here("analysis", "shelf_slope_gear_selectivity", "data", "catch_data_wide.rds")
)

model_list <- fit_list <- vector(mode = "list", length = length(unique_taxa))

for(ii in 1:length(unique_taxa)) {
  
  sel_spp <- unique_taxa[ii]
  
  message(ii, " ", sel_spp)
  
  fit_dat <- dplyr::filter(
    catch_data_wide, 
    COMMON_NAME == sel_spp,
    AREA_SWEPT_KM2_172>0, AREA_SWEPT_KM2_44>0
  )
  
  message(Sys.time(), ": Model 1")
  mod1 <- 
    glmmTMB::glmmTMB(
      formula = WEIGHT_172 ~ s(I(log(WEIGHT_44+1)), bs = "tp") + offset(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44)), 
      dispformula = ~ 1,
      data = fit_dat,
      family = glmmTMB::tweedie()
    )
  
  message(Sys.time(), ": Model 2")
  mod2 <- 
    glmmTMB::glmmTMB(
      formula = WEIGHT_172 ~ s(I(log(WEIGHT_44+1)), bs = "tp") + offset(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44)), 
      dispformula = ~ s(I(log(WEIGHT_44/AREA_SWEPT_KM2_44+1)), bs = "tp"),
      data = fit_dat,
      family = glmmTMB::lognormal()
    )
  
  message(Sys.time(), ": Model 3")
  mod3 <- 
    glmmTMB(
      formula = WEIGHT_172 ~ I(log(WEIGHT_44+1)) + offset(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44)), 
      dispformula = ~ 1,
      data = fit_dat,
      family = glmmTMB::tweedie()
    )
    
  message(Sys.time(), ": Model 4")
    mod4 <- 
      glmmTMB(
        formula = WEIGHT_172 ~ I(log(WEIGHT_44+1)) + offset(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44)), 
        dispformula = ~ s(I(log(WEIGHT_44/AREA_SWEPT_KM2_44+1)), bs = "tp"),
        data = fit_dat,
        family = glmmTMB::tweedie()
      )
    
    message(Sys.time(), ": Model 5")
    mod5 <- 
      glmmTMB(
        formula = WEIGHT_172 ~ 0 + I(log(WEIGHT_44+1)) + offset(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44)), 
        dispformula = ~ 1,
        data = fit_dat,
        family = glmmTMB::tweedie()
      )
    
    message(Sys.time(), ": Model 6")
    mod6 <- 
      glmmTMB(
        formula = WEIGHT_172 ~ 0 + I(log(WEIGHT_44+1)) + offset(log(AREA_SWEPT_KM2_172/AREA_SWEPT_KM2_44)), 
        dispformula = ~ s(I(log(WEIGHT_44/AREA_SWEPT_KM2_44+1)), bs = "tp"),
        data = fit_dat,
        family = glmmTMB::tweedie()
      )
    
    species_models <- list(mod1, mod2, mod3, mod4, mod5, mod6)
  
    aic_table <- make_aic_table(species_models)
    
    best_model <- species_models[[as.numeric(aic_table$model_name[1])]]
  
  model_list[[ii]] <- best_model 
  
  plot(DHARMa::simulateResiduals(model_list[[ii]]))
  
  fit_dat$fit <- 
    predict(
      model_list[[ii]],
      type = "link"
    )
  
  ex_fit <- 
    data.frame(
      WEIGHT_44 = exp(
        seq(
          log(ifelse(min(fit_dat$WEIGHT_44) == 0, min(fit_dat$WEIGHT_44)+0.0001, min(fit_dat$WEIGHT_44))),
          log(max(fit_dat$WEIGHT_44)),
          length = 300
        )),
      AREA_SWEPT_KM2_172 = 1,
      AREA_SWEPT_KM2_44 = 1,
      COMMON_NAME = sel_spp
    )
  
  ex_fit <- 
    cbind(
      ex_fit,
      as.data.frame(
        predict(model_list[[ii]], 
                newdata = ex_fit,
                type = "link",
                se.fit = TRUE
        )
      )
    )
  
  ex_fit$fit_response <- family(model_list[[ii]])$linkinv(ex_fit$fit)
  ex_fit$lwr_response <- family(model_list[[ii]])$linkinv(ex_fit$fit - 2 * ex_fit$se.fit)
  ex_fit$upr_response <- family(model_list[[ii]])$linkinv(ex_fit$fit + 2 * ex_fit$se.fit)
  
  fit_list[[ii]] <- ex_fit
  
}

fit_df <- do.call(what = rbind, args = fit_list)

p_gam_fits <- 
  ggplot() +
  geom_ribbon(
    data = fit_df,
    mapping = aes(x = WEIGHT_44, ymin = lwr_response, ymax = upr_response),
    alpha = 0.5
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_path(
    data = fit_df,
    mapping = aes(x = WEIGHT_44, y = fit_response)
  ) +
  geom_point(
    data = catch_data_wide,
    mapping = aes(x = WEIGHT_44,
                  y = WEIGHT_172)
  ) +
  scale_x_log10(name = "83-112 catch (kg)") +
  scale_y_log10(name = "PNE catch (kg)") +
  facet_wrap(~COMMON_NAME, scales = "free") +
  theme_bw()
