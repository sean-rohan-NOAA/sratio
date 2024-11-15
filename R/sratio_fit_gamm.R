#' Fit a logit-linked binomial or beta GAM with observation weights to catch comparison rate
#'
#' Fit a binomial or beta generalized additive model (GAM) to catch comparison rate (p) using `mgcv`.
#' Model can be fit with observations weighted by total count or model fits to residual variance, where residuals are modeled as a function of total count.
#'
#' @param data A data frame containing the data used for fitting the model. Must include catch comparison rate (p), size, and block (grouping variable for paired samples). A total_count field must be included when obs_weight_control$method is `"count"` or `"residuals_by_count"`.
#' @param k The basis dimension used for smooth terms in the GAM. 
#' @param gam_formula A formula for the GAM, by default `p ~ s(size, bs = "tp", k = k) + s(block, bs = "re")`.
#'   Smooth terms should be specified using `mgcv` syntax.
#' @param gam_family GAM family passed to `"mgcv::gam(family = gam_family)"` for the GAM. Defaults to `binomial(link = "logit")`.
#' @param obs_weight_control A list specifying control parameters for observation weighting. 
#'   Available options are:
#'   \describe{
#'     \item{`method`}{Character string specifying the weighting method. `"none"` (default) for no weights, `"count"` to weight by total count, and `"residuals_by_count"` to weight by model-estimated residual variance.}
#'     \item{`max_count`}{Maximum count allowed for each observation (numeric). Counts above this threshold are capped. Defaults = `Inf` (no threshold).}
#'     \item{`residual_type`}{Type of residuals to use for weighting model: `"absolute"` for absolute residuals or `"squared"` for squared residuals. Used only if `method` is `"residuals_by_count"`.}
#'     \item{`normalize_weight`}{Logical value indicating if weights should be normalized (divided by the mean weight). Normalizing weights allows for a consistent effective sample size and weighting scale among different weighting schemes. Default = `FALSE`.}
#'   }
#'
#' @return A fitted GAM object from `mgcv::gam`, containing the catch comparison rate model with the specified formula, family, and observation weights.
#'   
#' @import mgcv dplyr
#' @export

sratio_fit_gamm <- function(data,
                            k = 10,
                            gam_formula = 
                              p ~ s(size, bs = "tp", k = k) + s(block, bs = "re"),
                            gam_family = 
                              binomial(link = "logit"), 
                            obs_weight_control = 
                              list(method = "none",
                                   max_count = Inf,
                                   residual_type = NA,
                                   normalize_weight = FALSE)) {
  
  obs_weight <- 1
  resid_mod <- NA
  
  if(is.numeric(obs_weight_control$max_count)) {
    
    data$total_count <- 
      dplyr::if_else(
        data$total_count > obs_weight_control$max_count,
        obs_weight_control$max_count,
        data$total_count
      )
    
  }
  
  # Weight by raw counts
  if(obs_weight_control$method == "count") {
    
    obs_weight <- data$total_count
    
  }
  
  # Weight by model estimate of residual variance
  if(obs_weight_control$method == "residuals_by_count") { 
    
    mod <- 
      mgcv::gam(
        formula = gam_formula,
        data = data,
        family = gam_family
      )
    
    # Calculate absolute or squared residuals
    if(obs_weight_control$residual_type == "absolute") {
      data$transformed_resid <- abs(resid(mod))
    }
    
    if(obs_weight_control$residual_type == "squared") {
      data$transformed_resid <- resid(mod)^2
    }
    
    # Fit model model residual and predict
    resid_mod <- 
      mgcv::gam(
        formula = transformed_resid ~ s(total_count, k = 4, bs = "tp"), 
        data = data
      )
    
    # Estimate residual fit for each observation
    obs_weight <- 1/predict(object = resid_mod, newdata = data, type = "response")
    
    # Transform estimate of SD to variance when model is fit to absolute residuals  
    if(obs_weight_control$residual_type == "absolute") {
      
      obs_weight <- obs_weight^2
      
    }
    
  }
  
  # Normalize weights
  if(obs_weight_control$normalize_weight) {
    
    obs_weight <- obs_weight/mean(obs_weight)
  }
  
  data$obs_weight <- obs_weight
  
  mod <- 
    mgcv::gam(
      formula = gam_formula,
      data = data,
      weights = obs_weight,
      family = gam_family
    )
  
  return(list(mod = mod, 
              resid_mod = resid_mod))
  
}
