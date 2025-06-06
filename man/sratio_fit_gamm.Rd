% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sratio_fit_gamm.R
\name{sratio_fit_gamm}
\alias{sratio_fit_gamm}
\title{Fit a logit-linked binomial or beta GAM with observation weights to catch comparison rate}
\usage{
sratio_fit_gamm(
  data,
  k = 10,
  gam_formula = formula(p ~ s(size, bs = "tp", k = k) + s(block, bs = "re")),
  gam_family = binomial(link = "logit"),
  obs_weight_control = list(method = "none", max_count = Inf, residual_type = NA,
    normalize_weight = FALSE)
)
}
\arguments{
\item{data}{A data frame containing the data used for fitting the model. Must include catch comparison rate (p), size, and block (grouping variable for paired samples). A total_count field must be included when obs_weight_control$method is `"count"` or `"residuals_by_count"`.}

\item{k}{The basis dimension used for smooth terms in the GAM.}

\item{gam_formula}{A formula for the GAM, by default `p ~ s(size, bs = "tp", k = k) + s(block, bs = "re")`.
Smooth terms should be specified using `mgcv` syntax.}

\item{gam_family}{GAM family passed to `"mgcv::gam(family = gam_family)"` for the GAM. Defaults to `binomial(link = "logit")`.}

\item{obs_weight_control}{A list specifying control parameters for observation weighting. 
Available options are:
\describe{
  \item{`method`}{Character string specifying the weighting method. `"none"` (default) for no weights, `"count"` to weight by total count, and `"residuals_by_count"` to weight by model-estimated residual variance.}
  \item{`max_count`}{Maximum count allowed for each observation (numeric). Counts above this threshold are capped. Defaults = `Inf` (no threshold).}
  \item{`residual_type`}{Type of residuals to use for weighting model: `"absolute"` for absolute residuals or `"squared"` for squared residuals. Used only if `method` is `"residuals_by_count"`.}
  \item{`normalize_weight`}{Logical value indicating if weights should be normalized (divided by the mean weight). Normalizing weights allows for a consistent effective sample size and weighting scale among different weighting schemes. Default = `FALSE`.}
}}
}
\value{
A fitted GAM object from `mgcv::gam`, containing the catch comparison rate model with the specified formula, family, and observation weights.
}
\description{
Fit a binomial or beta generalized additive model (GAM) to catch comparison rate (p) using `mgcv`.
Model can be fit with observations weighted by total count or model fits to residual variance, where residuals are modeled as a function of total count.
}
