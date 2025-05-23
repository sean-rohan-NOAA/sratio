% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sratio_fit_bootstrap.R
\name{sratio_fit_bootstrap}
\alias{sratio_fit_bootstrap}
\title{Bootstrap selectivity ratio}
\usage{
sratio_fit_bootstrap(
  x,
  treatment_order,
  size_col,
  block_col,
  treatment_col,
  count_col,
  effort_col,
  sampling_factor_col = NULL,
  obs_weight_control = list(method = "none", max_count = Inf, residual_type = NA,
    normalize_weights = FALSE),
  gam_family,
  gam_formula = p ~ s(size, bs = "tp", k = k) + s(block, bs = "re"),
  k,
  scale_method = "sv",
  sratio_type = "absolute",
  n_cores = 1
)
}
\arguments{
\item{x}{list containing bootstrap samples, output from two_stage_bootstrap}

\item{treatment_order}{rank order of treatment levels, must match treatment names in treatment col}

\item{size_col}{Character vector; name of the size column}

\item{block_col}{Character vector; name of the sampling block column}

\item{treatment_col}{Character vector; name of the column containing treatment names}

\item{count_col}{Character vector; name of the column containing counts}

\item{effort_col}{Character vector; name of the column containing effort}

\item{sampling_factor_col}{Character vector; name of the column containing sampling factor (count multiplier, where estimated catch-at-length = count * sampling_factor)}

\item{obs_weight_control}{A list indicating which method to use to weight observations. See help documentation `sratio_fit_gamm()` for information about options (?sratio_fit_gamm) .}

\item{gam_family}{"binomial" or "beta"}

\item{gam_formula}{A formula for the GAM, by default `p ~ s(size, bs = "tp", k = k) + s(block, bs = "re")`.
Smooth terms should be specified using `mgcv` syntax.}

\item{k}{k for mgcv spline for size}

\item{scale_method}{Method to use for scaling the catch comparison rate for beta regression. See ?scale_for_betareg}

\item{sratio_type}{Which selectivity ratio calculation should be used? Absolute ("absolute") or relative ("relative")?}

\item{n_cores}{Number of cores to use for parallel processing.}
}
\description{
Uses the output of two_stage_bootstrap to fit binomial or beta generalized additive mixed models for catch comparison rate to generate selectivity ratio estimates.
}
